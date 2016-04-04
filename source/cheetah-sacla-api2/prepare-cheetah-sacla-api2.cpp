#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstdio>
#include <cmath>
#include <boost/format.hpp>
#include "DataAccessUserAPI.h"
#include "hdf5.h"

static std::string get_geom(int run);
static void get_geom_h5(int run);
static int add_image(double *buffer, int tag, double photon_energy);
static int run(int runid);

const int bl = 3;
const int xsize = 512;
const int ysize = 1024; 
const int ydatasize = 1030;
const int blocksize = xsize * ysize;
const int buffersize = blocksize * 8;
const int stride = 2; // 30 Hz mode (60 / 30)
char *det_name_template[30] = {"EventInfo_stor0%d", "MPCCD-8-2-001-%d", "EventInfo_stor1_0%d", "MPCCD-8-2-002-%d"};
char *LLF_ID = "BL3-ST4-MPCCD-octal"; // FIXME: LLF ID changed to ST4. how to switch automatically??

// FIXME: make these local
double buffer[buffersize] = {};
unsigned short averaged[buffersize] = {};
float posx[buffersize] = {}, posy[buffersize] = {}, posz[buffersize] = {};
float gains[9] = {};
int det_temp_idx = -1;
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

int main(int argc, char **argv) {
  if(argc < 2) {
    std::cout << argv[0] << " RunID" << std::endl;
    return 1;
  }
  
  int runid = atoi(argv[1]);
  run(runid);
  
  return 0;
}

int run(int runid) {
  int retno = 0, start, end, tag_hi = -1;

  printf("SACLA geometry & dark average exporter\n");
  printf(" By Takanori Nakane\n");
  printf(" version 20160127 with new API\n\n");
  
  // Get tag_hi, start, end
  retno = sy_read_start_tagnumber(&tag_hi, &start, bl, runid);
  retno += sy_read_end_tagnumber(&tag_hi, &end, bl, runid);
  if (retno != 0) {
	printf("ERROR: Cannot read run %d.\n", runid);
	printf("If this run is before Nov 2014, please use the old API version.\n");
    return -1;
  }

  int numAll = (end - start) / stride + 1;
  printf("Run %d contains tag %d - %d (%d images)\n", runid, start, end, numAll);
  int *tagAll = (int*)malloc(sizeof(int) * numAll);
  for (int i = 0; i < numAll; i++) {
    tagAll[i] = start + i * stride;
  }

  // How many dark frames?
  std::vector<std::string> shutter;
  if (myReadSyncDataList(&shutter, "xfel_bl_3_shutter_1_open_valid/status", tag_hi, numAll, tagAll) != 0) {
    printf("Failed to get shutter status.\n");
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

  int *tagList = (int*)malloc(sizeof(int) * numDark);
  for (int i = 0; i < numDark; i++) {
    tagList[i] = start + i * stride;
  }

  // find detector ID
  struct da_string_array *det_ids;
  int n_detid;

  printf("Detector configulation:\n");
  da_alloc_string_array(&det_ids);
  sy_read_detidlist(det_ids, bl, runid);
  
  da_getsize_string_array(&n_detid, det_ids);
  for (int i = 0; i < n_detid; i++) {
    char *detid;
    da_getstring_string_array(&detid, det_ids, i);
    printf(" detID #%d = %s\n", i, detid);
    if (strcmp(detid, "MPCCD-8-2-001-1") == 0) {
      det_temp_idx = 1;
    } else if (strcmp(detid, "MPCCD-8-2-002-1") == 0) {
      det_temp_idx = 3;
    } else if (strcmp(detid, "EventInfo_stor01") == 0) {
      printf("ERROR: This detector is not longer supported by the API. Use old Cheetah.\n");
	}
    free(detid);
  }
  if (det_temp_idx == -1) {
    printf("ERROR: Unknown or non-supported detector ID.\n");
    return -1;
  }
  da_destroy_string_array(&det_ids);

  printf("\n");
   
  // Create storage readers and buffers
  // collect first image (used to generate geometry files)
  printf("Initializing reader and buffer\n");
  for (int det_id = 0; det_id < 8; det_id++) {
    char det_name[256];
    snprintf(det_name, 256, det_name_template[det_temp_idx], det_id + 1);

    printf(" detector %s\n", det_name);
    retno = st_create_streader(&streaders[det_id], det_name, bl, 1, &runid);
    if (retno != 0) {
      printf("Failed to create streader for %s.\n", det_name);
      return -1;
    }
    retno = st_create_stbuf(&databufs[det_id], streaders[det_id]);
    if (retno != 0) {
      printf("Failed to allocate databuffer for %s.\n", det_name);
      return -1;
    }
    uint32_t tagid = start;
    retno = st_collect_data(databufs[det_id], streaders[det_id], &tagid);
    if (retno != 0) {
      printf("Failed to collect data for %s.\n", det_name);
      return -1;
    }
  }

  // Geometry generation
  printf("\n");
  char filename[256];
  snprintf(filename, 256, "%06d.geom", runid);
  std::ofstream ofs(filename);
  ofs << get_geom(runid);
  ofs.close();
  printf("CrystFEL geometry file was written to %s\n", filename);
  get_geom_h5(runid);

  // LLF values
  printf("\n");
  struct da_string_array *llf;
  int n_llf;
  
  da_alloc_string_array(&llf);
  sy_read_statistics_detllf(llf, LLF_ID, bl, tag_hi, numDark, tagList);

  da_getsize_string_array(&n_llf, llf);
  printf("\nLLF statistics for %d frames.\n", n_llf);
  for (int i = 0; i < n_llf; i++) {
    char *val;
    da_getstring_string_array(&val, llf, i);
    printf(" %s", val);
    free(val);
  }
  printf("\n\n");
  da_destroy_string_array(&llf);

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
  std::vector<std::string> pulse_energies;
  struct da_string_array *energies;
  int n_energy;
  
  if (runid >=333661 && runid <= 333682) {
    // 19 May 2015: spectrometer broken! use config value instead
    for (int i = 0; i < numDark; i++) {
      pulse_energies.push_back("7.0");
    }
  } else {
    da_alloc_string_array(&energies);
    sy_read_syncdatalist(energies, "xfel_bl_3_tc_spec_1/energy", tag_hi, numDark, tagList);

    da_getsize_string_array(&n_energy, energies);
    for (int i = 0; i < n_energy; i++) {
      char *val;
      da_getstring_string_array(&val, energies, i);
      pulse_energies.push_back(val);
      free(val);
    }
    da_destroy_string_array(&energies);
  }
  
  printf("xfel_bl_3_tc_spec_1/energy for %ld frames.\n", pulse_energies.size());
  for (unsigned int i = 0; i < pulse_energies.size(); i++) {
    printf(" %s", pulse_energies[i].c_str());
  }
  printf("\n\n");

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
    num_added += add_image(buffer, tagID, 1000 * atof(pulse_energies[j].c_str()));
  }
  printf("\nDone. Averaged %d frames.\n", num_added);
  
  if (num_added < 1) return -1;
  for (int i = 0; i < buffersize; i++) {
    double tmp = buffer[i] / num_added;
    if (tmp < 0) averaged[i] = 0;
    if (tmp > USHRT_MAX) averaged[i] = USHRT_MAX;
    else averaged[i] = (unsigned short)tmp;
    // TODO: Is this correct treatmen2Dt?
  }
  
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
  
  hid_t file_id, dataset_id, dataspace_id, group_id;
  hsize_t dims[] = {ysize * 8, xsize};
  snprintf(filename, 256, "%06d-dark.h5", runid);
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple(2, dims, NULL);
  group_id = H5Gcreate2(file_id, "/data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataset_id = H5Dcreate2(file_id, "/data/data", H5T_NATIVE_USHORT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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

static std::string get_geom(int run) {
  std::stringstream ss;
  float detx, dety, detz, rotation;
//  float stage_dir, shift_weight, aperture_origin, aperture_par;
//  int manipulator_pos;

  double energy;
  sy_read_config_photonenergy(&energy, bl, run); // returned in keV

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
     << boost::format("photon_energy = /%%/photon_energy_ev ; roughly %f. change this to /LCLS/photon_energy_eV for single-event files\n") % (energy * 1000)
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
    ss << boost::format("q%d/adu_per_eV = %f\n") % panel % (1.0 / (0.1 * energy * 1000)); // Keitaro's 0.1 photon
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
  hid_t file_id, dataset_id, dataspace_id;
  hsize_t dims[] = {ysize * 8, xsize};

  snprintf(filename, 256, "%d-geom.h5", runid);
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple(2, dims, NULL);
  dataset_id = H5Dcreate2(file_id, "x", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posx);
  dataset_id = H5Dcreate2(file_id, "y", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posy);
  dataset_id = H5Dcreate2(file_id, "z", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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

