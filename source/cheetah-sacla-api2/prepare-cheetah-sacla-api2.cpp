#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <boost/format.hpp>
#include "DataAccessUserAPI.h"
#include "hdf5.h"

static std::string get_geom(int run, double energy);
static void get_geom_h5(int run);
static int add_image(double *buffer, int tag, double photon_energy);
static int run(int runid);

int bl = 3;
const int xsize = 512;
const int ysize = 1024; 
const int ydatasize = 1030;
const int blocksize = xsize * ysize;
const int buffersize = blocksize * 8;
const char *LLF_ID = "BL3-ST4-MPCCD-octal"; // FIXME: LLF ID changed to ST4. how to switch automatically??

// FIXME: make these local
double buffer[buffersize] = {};
unsigned short averaged[buffersize] = {};
float posx[buffersize] = {}, posy[buffersize] = {}, posz[buffersize] = {};
float gains[9] = {};
char *streaders[8], *databufs[8];

int myReadSyncDataList(std::vector<std::string>* buffer, char* equip_id, int tag_hi, int n_taglist, int *taglist) {
	struct da_string_array *strbuf;
	int n_strbuf, retno;

	da_alloc_string_array(&strbuf);
	retno = sy_read_syncdatalist(strbuf, equip_id, tag_hi, n_taglist, taglist);
	if (retno != 0) {
		return retno;
	}

	da_getsize_string_array(&n_strbuf, strbuf);
	for (int i = 0; i < n_strbuf; i++) {
		char *val;
		da_getstring_string_array(&val, strbuf, i);
		buffer->push_back(val);
		free(val);
	}
	da_destroy_string_array(&strbuf);

	return retno;
}

void log_error(char *message) {
  FILE *status = fopen("status.txt", "w");
  fprintf(status, "Status: Status=Error-%s\n", message);
  fclose(status);
}

/*
 * Compression helps a lot!
 *  2016-Dec BL2 Run 21498
 *
 *  dark average 8391040, deflate level4 216947, shuf + deflate 179189
 *  geometry    50333792,              20818705,               1490293
 *
 */

int main(int argc, char **argv) {
  if(argc < 2) {
    std::cout << argv[0] << " RunID [Beamline (2 or 3)]" << std::endl;
    return -1;
  }
  
  int runid = atoi(argv[1]);
  if (argc >= 3) {
    bl = atoi(argv[2]);
    if (bl != 2 && bl != 3) {
      printf("Beamline must be 2 or 3.\n");
      return -1;
    }
  }
  printf("Set beamline to %d\n", bl);
  run(runid);
  
  return 0;
}

int run(int runid) {
  int retno = 0, start, end, tag_hi = -1;

  printf("SACLA geometry & dark average exporter\n");
  printf(" By Takanori Nakane\n");
  printf(" version 20170209 with new API\n\n");
  printf("WARNING: This C++ version is no longer maintained. Use prepare-cheetah-sacla-api2.py!\n\n");
  
  // Get tag_hi, start, end
  retno = sy_read_start_tagnumber(&tag_hi, &start, bl, runid);
  retno += sy_read_end_tagnumber(&tag_hi, &end, bl, runid);
  if (retno != 0) {
    printf("ERROR: Cannot read run %d.\n", runid);
    printf("If this run is before Nov 2014, please use the old API version.\n");
    log_error("BadRunID");
    return -1;
  }

  struct da_int_array *tagbuf;
  int numAll;
  da_alloc_int_array(&tagbuf, 0, NULL);
  sy_read_taglist_byrun(tagbuf, bl, runid);
  da_getsize_int_array(&numAll, tagbuf);
  printf("Run %d contains tag %d (inclusive) to %d (exclusive), thus %d images\n", runid, start, end, numAll);
  int *tagAll = (int*)malloc(sizeof(int) * numAll);
  for (int i = 0; i < numAll; i++) {
    da_getint_int_array(tagAll + i, tagbuf, i);
  }

  // How many dark frames?
  std::vector<std::string> shutter;
  if ((bl == 3 && myReadSyncDataList(&shutter, "xfel_bl_3_shutter_1_open_valid/status", tag_hi, numAll, tagAll) != 0) ||
      (bl == 2 && myReadSyncDataList(&shutter, "xfel_bl_2_shutter_1_open_valid/status", tag_hi, numAll, tagAll) != 0)) {
    printf("Failed to get shutter status.\n");
    log_error("NoShutterStatus");
    return -1;
  }
  free(tagAll);

  int numDark = 0;
  for (int i = 0; i < numAll; i++) {
    if (atoi(shutter[i].c_str()) == 0) {
      numDark++;
    } else {
      break;
    }
  }
  printf("Number of dark frames: %d\n\n", numDark);
  if (numDark == 0) {
    printf("ERROR: No dark image!\n");
    log_error("NoDarkImage");
  }

  int *tagList = (int*)malloc(sizeof(int) * numDark);
  for (int i = 0; i < numDark; i++) {
    da_getint_int_array(tagList + i, tagbuf, i);
  }
  da_destroy_int_array(&tagbuf);

  // find detector ID
  struct da_string_array *det_ids;
  int n_detid;
  char det_template[256] = {};

  printf("Detector configulation:\n");
  da_alloc_string_array(&det_ids);
  sy_read_detidlist(det_ids, bl, runid);
  
  da_getsize_string_array(&n_detid, det_ids);
  for (int i = 0; i < n_detid; i++) {
    char *detid;
    da_getstring_string_array(&detid, det_ids, i);
    printf(" detID #%02d = %s\n", i, detid);
    if (strncmp(detid, "MPCCD-8", 7) == 0) {
      int len = strlen(detid);
      if (detid[len - 2] == '-' &&  detid[len - 1] == '1') {
        printf("  prefix and suffix matched. using this as the detector name template.\n");
        strncpy(det_template, detid, 255);
      } 
    }
    if (strcmp(detid, "EventInfo_stor01") == 0) {
      printf("ERROR: This detector is not longer supported by the API. Use old Cheetah.\n");
    }
    free(detid);
  }
  if (det_template[0] == 0) {
    printf("ERROR: Unknown or non-supported detector ID.\n");
    log_error("NoSupportedDetectorFound");
    return -1;
  }
  da_destroy_string_array(&det_ids);

  printf("\n");
   
  // Create storage readers and buffers
  // collect first image (used to generate geometry files)
  printf("Initializing reader and buffer:\n");
  for (int det_id = 0; det_id < 8; det_id++) {
    char det_name[256];
    strncpy(det_name, det_template, 255);
    det_name[strlen(det_name) - 1] = '0' + det_id + 1;

    printf(" detector %s\n", det_name);
    retno = st_create_streader(&streaders[det_id], det_name, bl, 1, &runid);
    if (retno != 0) {
      printf("Failed to create streader for %s.\n", det_name);
      log_error("FailedOn_create_streader");
      return -1;
    }
    retno = st_create_stbuf(&databufs[det_id], streaders[det_id]);
    if (retno != 0) {
      printf("Failed to allocate databuffer for %s.\n", det_name);
      log_error("FailedOn_create_stbuf");
      return -1;
    }
    uint32_t tagid = start;
    retno = st_collect_data(databufs[det_id], streaders[det_id], &tagid);
    if (retno != 0) {
      printf("Failed to collect data for %s.\n", det_name);
      log_error("FailedOn_collect_data");
      return -1;
    }
  }

  // PD values
  struct da_string_array *pd_laser;
  int n_pd_laser;
  const char *laser_id = "xfel_bl_3_st_4_pd_laser_fitting_peak/voltage";
  
  da_alloc_string_array(&pd_laser);
  sy_read_syncdatalist(pd_laser, laser_id, tag_hi, numDark, tagList);

  da_getsize_string_array(&n_pd_laser, pd_laser);
  printf("%s for %d frames.\n", laser_id, n_pd_laser);
  for (int i = 0; i < n_pd_laser; i++) {
    char *val;
    da_getstring_string_array(&val, pd_laser, i);
    printf(" %s", val);
    free(val);
  }
  printf("\n\n");
  da_destroy_string_array(&pd_laser);

  // Pulse energies
  std::vector<float> pulse_energies;
  struct da_string_array *energies;
  int n_energy = 0;
  float mean_energy = 0;
  
  if (runid >=333661 && runid <= 333682) {
    // 19 May 2015: spectrometer broken! use config value instead
    for (int i = 0; i < numDark; i++) {
      pulse_energies.push_back(7.0);
    }
  } else {
    da_alloc_string_array(&energies);
    if (bl == 3) {
      sy_read_syncdatalist(energies, "xfel_bl_3_tc_spec_1/energy", tag_hi, numDark, tagList);
    } else if (bl == 2) {
      sy_read_syncdatalist(energies, "xfel_bl_2_tc_spec_1/energy", tag_hi, numDark, tagList);
    }

    da_getsize_string_array(&n_energy, energies);
    for (int i = 0; i < n_energy; i++) {
      char *val;
      da_getstring_string_array(&val, energies, i);
      pulse_energies.push_back(1000.0 * atof(val));
      free(val);
    }
    da_destroy_string_array(&energies);
  }
  
  printf("Photon energies (eV) for %d frames:\n", n_energy);
  for (unsigned int i = 0; i < n_energy; i++) {
    printf(" %.2f", pulse_energies[i]);
    mean_energy += pulse_energies[i];
  }
  mean_energy /= n_energy;
  printf("\nMean photon energy = %.2f eV\n", mean_energy);
  double config_energy;
  sy_read_config_photonenergy(&config_energy, bl, runid);
  printf("Accelerator config = %.2f eV\n\n", 1000.0 * config_energy);

  // Geometry generation
  char filename[256];
  snprintf(filename, 256, "%d.geom", runid);
  std::ofstream ofs(filename);
  ofs << get_geom(runid, mean_energy);
  ofs.close();
  printf("CrystFEL geometry file was written to %s\n", filename);
  get_geom_h5(runid);
  printf("\n");

  // Dark averages
  
  int num_added = 0;
  for (int j = 0; j < numDark; j++) {
    int tagID = tagList[j];
    printf("Processing tag %d (%.1f%% done)\n", tagID, 100.0 * (j + 1) / numDark);
    if (j % 3 == 0) {
      FILE *status = fopen("status.txt", "w");
      fprintf(status, "Status: Total=%d,Processed=%d,Status=DarkAveraging\n", numDark, j + 1);
      fclose(status);
    }
    num_added += add_image(buffer, tagID, pulse_energies[j]);
  }
  printf("\nDone. Averaged %d frames.\n", num_added);
 
  if (num_added < 1) {
    log_error("NoImageAveraged");
    return -1;
  }
  int n_neg = 0, n_overflow = 0;
  for (int i = 0; i < buffersize; i++) {
    double tmp = buffer[i] / num_added;
    if (tmp < 0) {averaged[i] = 0; n_neg++;}
    else if (tmp > USHRT_MAX) {averaged[i] = USHRT_MAX; n_overflow++;}
    else averaged[i] = (unsigned short)tmp;
  }
  printf(" #neg (< 0) %d #overflow (> %d) %d\n", n_neg, USHRT_MAX, n_overflow);
  
  /*
    snprintf(filename, 256, "%d-dark.img", runid);
    FILE *fh = fopen(filename, "wb");
    fprintf(fh, "{\n"
    "HEADER_BYTES=512;\n"
    "DIM=2;\n"
    "BYTE_ORDER=little_endian;\n"
    "TYPE=unsigned_short;\n"
    "SIZE1=%d;\n"
    "SIZE2=%d;\n"
    "PIXEL_SIZE=%f;\n"
    "WAVELENGTH=%f;\n"
    "DISTANCE=%f;\n"
    "PHI=0.0;\n"
    "OSC_START=0.00;\n"
    "OSC_END=0.00;\n"
    "OSC_RANGE=0.00;\n"
    "AXIS=phi;\n"
    "BEAM_CENTER_X=%f;\n"
    "BEAM_CENTER_Y=%f;\n"
    "}\n",
    xsize, ysize * 8, 0.05,
    1.7, 515, 0, 0);
    fseek(fh, 512, SEEK_SET);
    fwrite(averaged, sizeof(unsigned short), blocksize * 8, fh);
    fclose(fh);
  */
  
  hid_t file_id, dataset_id, dataspace_id, group_id, plist_id;
  hsize_t dims[] = {ysize * 8, xsize};
  snprintf(filename, 256, "%d-dark.h5", runid);
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple(2, dims, NULL);
  group_id = H5Gcreate2(file_id, "/data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_shuffle(plist_id);
  H5Pset_deflate(plist_id, 4);
  H5Pset_chunk(plist_id, 2, dims);
  dataset_id = H5Dcreate2(file_id, "/data/data", H5T_NATIVE_USHORT, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, averaged);
  H5Sclose(dataspace_id);
  H5Gclose(group_id);
  H5Dclose(dataset_id);
  H5Fclose(file_id);
  printf("Dark average was written to %s\n", filename);

  for (int det_id = 0; det_id < 8; det_id++) {
    st_destroy_stbuf(&databufs[det_id]);
    st_destroy_streader(&streaders[det_id]);
  }

  free(tagList);
  return 0;
}

static std::string get_geom(int run, double energy) {
  std::stringstream ss;
  float detx, dety, detz, rotation;
//  float stage_dir, shift_weight, aperture_origin, aperture_par;
//  int manipulator_pos;

  ss << "; CrystFEL geometry file produced by prepare-cheetah-sacla-api2\n"
     << ";   by Takanori Nakane\n"
     << "; for tiled but NOT reassembled images (512x8192 pixels)\n\n"
     << "; NOTE:\n"
     << "; This file is for multi-event HDF5 files. To use in single-event\n"
     << "; files, you have to edit 'data' record and prepare .beam file\n\n"
     << "clen = 0.0515  ; 51.5 mm camera length. You SHOULD optimize this!\n"
     << "res = 20000  ; 50 micron  1m /50 micron\n"
     << ";badrow_direction = x\n"
     << ";adu_per_eV = /LCLS/adu_per_eV\n"
     << ";max_adu = 250000 ; should NOT be used. see best practice on Web\n"
     << "data = /%/data ; change this to /LCLS/data for single-event files\n"
     << boost::format("photon_energy = /%%/photon_energy_ev ; roughly %f. change this to /LCLS/photon_energy_eV for single-event files\n") % energy
     << "rigid_group_q1 = q1\n" 
     << "rigid_group_q2 = q2\n" 
     << "rigid_group_q3 = q3\n" 
     << "rigid_group_q4 = q4\n" 
     << "rigid_group_q5 = q5\n" 
     << "rigid_group_q6 = q6\n" 
     << "rigid_group_q7 = q7\n" 
     << "rigid_group_q8 = q8\n\n" 
     << "rigid_group_collection_connected = q1,q2,q3,q4,q5,q6,q7,q8\n"
     << "rigid_group_collection_independent = q1,q2,q3,q4,q5,q6,q7,q8\n\n";

  printf("Panel geometry:\n");

  for (int det_id = 0; det_id < 8; det_id++) {
    mp_read_absgain(&gains[det_id], databufs[det_id]);
    mp_read_posx(&detx, databufs[det_id]);
    mp_read_posy(&dety, databufs[det_id]);
    mp_read_posz(&detz, databufs[det_id]);
    mp_read_rotationangle(&rotation, databufs[det_id]);
    
    printf(" panel %d gain %f pos (%f, %f, %f) rotation %f energy %f\n", det_id, gains[det_id], detx, dety, detz, rotation, energy);

    /*
    mp_read_stagedirection(&stage_dir, databufs[det_id]);
    mp_read_shiftweight(&shift_weight, databufs[det_id]);
    ReadDetApertureOrigin(aperture_origin, det_name, bl, taghi, tag); // ?
    ReadDetAperturePar(aperture_par, det_name, bl, taghi, tag); // ?
    mp_read_manipulator_shift(&manipulator_pos, databufs[det_id]);
    //	  printf(" stage_dir %f shift_weight %f aperture_origin %f aperture_par %f manipulator_pos %d\n",
    //           stage_dir, shift_weight, aperture_origin, aperture_par, manipulator_pos);
    */
  
    // Thanks to Keitaro Yamashita
    detx /= 50; dety /= 50;
    rotation *= M_PI / 180;
    
    int row = det_id, panel = det_id + 1;

    // Nphotons = S [ADU] * G [e-/ADU] / (E [eV] * 3.65 [eV/e-]) according to the manual.
    // Thus, ADU/eV = 1/(3.65*G)

    //    ss << boost::format("q%d/adu_per_eV = %f\n") % panel % (1.0 / (gain * 3.65));
    ss << boost::format("q%d/adu_per_eV = %f\n") % panel % (1.0 / (0.1 * energy)); // Keitaro's 0.1 photon
    ss << boost::format("q%d/min_fs = %d\n") % panel % 0;
    ss << boost::format("q%d/min_ss = %d\n") % panel % (row * ysize);
    ss << boost::format("q%d/max_fs = %d\n") % panel % 511;
    ss << boost::format("q%d/max_ss = %d\n") % panel % ((row + 1) * ysize - 1);
    
    ss << boost::format("q%d/fs = %fx %+fy\n") % panel % std::cos(rotation) % std::sin(rotation);
    ss << boost::format("q%d/ss = %fx %+fy\n") % panel % -std::sin(rotation) % std::cos(rotation);
    ss << boost::format("q%d/corner_x = %f\n") % panel % -detx;
    ss << boost::format("q%d/corner_y = %f\n\n") % panel % dety;
  }
  return ss.str();
}

static void get_geom_h5(int runid) {
  float gain, detx, dety, detz, rotation;

  for (int det_id = 0; det_id < 8; det_id++) {
    mp_read_absgain(&gain, databufs[det_id]);
    mp_read_posx(&detx, databufs[det_id]);
    mp_read_posy(&dety, databufs[det_id]);
    mp_read_posz(&detz, databufs[det_id]);
    mp_read_rotationangle(&rotation, databufs[det_id]);

    // Thanks to Keitaro Yamashita
    rotation *= M_PI / 180;
    
    int row = det_id;
    double pixel_size = 50 * 1E-6; // 50um in m
    double fast_x = std::cos(rotation) * pixel_size, fast_y = std::sin(rotation) * pixel_size;
    double slow_x = -std::sin(rotation) * pixel_size, slow_y = std::cos(rotation) * pixel_size;

    double orgx = detx * 1E-6, orgy = dety * 1E-6; // um -> m
    
    for (int x = 0; x < xsize; x++) {
      for (int y = 0; y < ysize; y++) {
	int bufx = x, bufy = row * ysize + y;
	posx[bufx + xsize * bufy] = fast_x * x + slow_x * y - orgx;
	posy[bufx + xsize * bufy] = -(fast_y * x + slow_y * y + orgy);
	posz[bufx + xsize * bufy] = 0; // offset
      }
    }
  }

  char filename[256];
  hid_t file_id, dataset_id, dataspace_id, plist_id;
  hsize_t dims[] = {ysize * 8, xsize};

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_shuffle(plist_id);
  H5Pset_deflate(plist_id, 4);
  H5Pset_chunk(plist_id, 2, dims);

  snprintf(filename, 256, "%d-geom.h5", runid);
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple(2, dims, NULL);
  dataset_id = H5Dcreate2(file_id, "x", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posx);
  dataset_id = H5Dcreate2(file_id, "y", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posy);
  dataset_id = H5Dcreate2(file_id, "z", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posz);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);
  H5Fclose(file_id);
  printf("Cheetah geometry file was written to %s\n", filename);
}

 static int add_image(double *buffer, int tag, double photon_energy) {
  int retno = 0;
  float buf_panel[xsize * ydatasize];
  float gain;
  
  for (int det_id = 0; det_id < 8; det_id++) {
    uint32_t tagid = tag; // tagid will be overwritten
    retno = st_collect_data(databufs[det_id], streaders[det_id], &tagid);
    if (retno != 0 || tagid != (uint32_t)tag) {
      printf("WARNING: Failed to collect tag %d for detector #%det_id (0-indexed)\n", tag, det_id);
      return 0;
    }

    retno = st_read_det_data(buf_panel, databufs[det_id], 0);
    if (retno != 0) {
      printf("WARNING: Failed to read tag %d for detector #%det_id (0-indexed)\n", tag, det_id);
      return 0;
    }

    gain = gains[det_id] / gains[0];
    gain *= gains[0] * 3.65 / 0.1 / photon_energy; // CHECKME: Keitaro's 0.1 photon unit

    // Mark origin for debugging
    if (false) {
      for (int i = 0; i < 10; i++) {
	for (int j = 0; j < 10; j++) {
	  buf_panel[i * xsize + j] = det_id * 500;
	}
      }
    }
    
    int offset = det_id * blocksize;
    for (int i = 0; i < blocksize; i++) {
      buffer[offset + i] += buf_panel[i] * gain;
    }
  }
  
  return 1;
}

