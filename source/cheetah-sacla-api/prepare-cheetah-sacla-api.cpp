#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstdio>
#include <cmath>
#include <boost/format.hpp>
#include "SaclaDataAccessUserAPI.h"
#include "hdf5.h"

static std::string get_geom(int run, int taghi, int tag);
static void get_geom_h5(int runid, int taghi, int tag);
static int add_image(double *buffer, int run, int taghi, int tag, double photon_energy);
static int run(int runid);

const int bl = 3;
const int xsize = 512;
const int ysize = 1024; 
const int ydatasize = 1030;
const int blocksize = xsize * ysize;
const int buffersize = blocksize * 8;
const int numDark = 150;
const int stride = 2; // 30Hz mode (60 / 30)
char *det_name_template[30] = {"EventInfo_stor0%d", "MPCCD-8-2-001-%d", "EventInfo_stor1_0%d"};
char *LLF_ID = "BL3-ST4-MPCCD-octal"; // FIXME: LLF ID changed to ST4. how to switch automatically??

// FIXME: make these local
double buffer[buffersize] = {};
unsigned short averaged[buffersize] = {};
float posx[buffersize] = {}, posy[buffersize] = {}, posz[buffersize] = {};
float gains[9] = {};
int det_temp_idx = -1;

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
  int retno = 0, start;
  int tag_hi = -1;

  printf("SACLA geometry & dark average exporter\n");
  printf(" By Takanori Nakane \n\n");
  printf(" version 20151002 with old API\n\n");
 
  // get tag_hi and start
  ReadStartTagNumber(tag_hi, start, bl, runid);
  std::vector<std::string> det_ids;
  ReadDetIDList(&det_ids, bl, runid);  

  for (unsigned int i = 0; i < det_ids.size(); i++) {
    if (det_ids[i] == "MPCCD-8-2-001-1") {
      det_temp_idx = 1;
      break;
    } else if (det_ids[i] == "EventInfo_stor01") {
      det_temp_idx = 0;
      break;
    } else if (det_ids[i] == "EventInfo_stor1_01") { // FIXME: not working?
      det_temp_idx = 2;
      break;
    }
  }
  if (det_temp_idx == -1) {
    printf("ERROR: Unknown detector ID %s.\n", det_ids[0].c_str());
    return -1;
  }

  std::vector<int> tagList;
  for (int i = 0; i < numDark; i++)
    tagList.push_back(start + i * stride);
  //runid = 209050; tagList.clear(); tagList.push_back(121943650); // for debugging
  
  char filename[256];
  snprintf(filename, 256, "%06d.geom", runid);
  std::ofstream ofs(filename);
  ofs << get_geom(runid, tag_hi, tagList[0]);
  ofs.close();
  get_geom_h5(runid, tag_hi, tagList[0]);
  printf("CrystFEL geometry file was written to %s\n", filename);

  // TODO: is this right??
  std::vector<std::string> llf;
  ReadStatisticsOfDetLLF(&llf, LLF_ID, bl, tag_hi, tagList);
  printf("\nLLF statistics for %ld frames.\n", llf.size());
  for (unsigned int i = 0; i < llf.size(); i++)
    printf(" %s", llf[i].c_str());
  printf("\n");

  std::vector<std::string> pd_laser, pd_user2;
  retno = ReadSyncDataList(&pd_laser, "xfel_bl_3_st_4_pd_laser_fitting_peak/voltage", tag_hi, tagList);
  printf("\nxfel_bl_3_st_4_pd_laser_fitting_peak/voltage for %ld frames.\n", pd_laser.size());
  for (unsigned int i = 0; i < pd_laser.size(); i++)
    printf(" %s", pd_laser[i].c_str());
  printf("\n\n");
  retno = ReadSyncDataList(&pd_user2, "xfel_bl_3_st_4_pd_user_2_fitting_peak/voltage", tag_hi, tagList);
  printf("\nxfel_bl_3_st_4_pd_user_2_fitting_peak/voltage for %ld frames.\n", pd_user2.size());
  for (unsigned int i = 0; i < pd_user2.size(); i++)
    printf(" %s", pd_user2[i].c_str());
  printf("\n\n");

  std::vector<std::string> pulse_energies;
  if (runid >=333661 && runid <= 333682) {
    // 19 May 2015: spectrometer broken! use config value instead
    for (unsigned int i = 0; i < tagList.size(); i++) {
      pulse_energies.push_back("7.0");
    }
  } else {
    retno = ReadSyncDataList(&pulse_energies, "xfel_bl_3_tc_spec_1/energy", tag_hi, tagList);
  }
  printf("\nxfel_bl_3_tc_spec_1/energy for %ld frames.\n", pulse_energies.size());

  for (unsigned int i = 0; i < pulse_energies.size(); i++)
    printf(" %s", pulse_energies[i].c_str());
  printf("\n\n");
  
  int num_added = 0;
  for (unsigned int j = 0, tagSize = tagList.size(); j < tagSize; j++) {
    int tagID = tagList[j];
    printf("Processing tag %d (%.1f%% done)\n", tagID, 100.0 * (j + 1) / tagSize);
    if (j % 3 == 0) {
      FILE *status = fopen("status.txt", "w");
      fprintf(status, "Status: Total=%d,Processed=%d,Status=DarkAveraging\n", tagSize, j + 1);
      fclose(status);
    }
    num_added += add_image(buffer, runid, tag_hi, tagID, 1000 * atof(pulse_energies[j].c_str()));
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
  return 0;
}

static std::string get_geom(int run, int taghi, int tag) {
	std::stringstream ss;
	float gain, detx, dety, detz, rotation;
//  float stage_dir, shift_weight, aperture_origin, aperture_par;
//  int manipulator_pos;

  double energy;
  ReadConfigOfPhotonEnergy(energy, bl, taghi, tag); // returned in keV

  ss << "; CrystFEL geometry file produced by averageImageGain\n"
     << ";   Takanori Nakane (takanori.nakane@bs.s.u-tokyo.ac.jp)\n"
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
  
  char det_name[256];
  for (int det_id = 1; det_id <= 8; det_id++) {
    snprintf(det_name, 256, det_name_template[det_temp_idx], det_id);
    
    ReadAbsGain(gain, det_name, bl, run, taghi, tag);
    ReadDetPosX(detx, det_name, bl, run, taghi, tag);
    ReadDetPosY(dety, det_name, bl, run, taghi, tag);
    ReadDetPosZ(detz, det_name, bl, run, taghi, tag);
    ReadDetRotationAngle(rotation, det_name, bl, run, taghi, tag);
    
    //    printf("gain %f pos (%f, %f, %f) rotation %f energy %f\n", gain, detx, dety, detz, rotation, energy);
    /*
    ReadDetStageDirection(stage_dir, det_name, bl, taghi, tag);
    ReadDetShiftWeight(shift_weight, det_name, bl, taghi, tag);
    ReadDetApertureOrigin(aperture_origin, det_name, bl, taghi, tag);
    ReadDetAperturePar(aperture_par, det_name, bl, taghi, tag);
    ReadDetManipulatorPos(manipulator_pos, det_name, bl, taghi, tag);
    //	  printf("stage_dir %f shift_weight %f aperture_origin %f aperture_par %f manipulator_pos %d\n",
    //           stage_dir, shift_weight, aperture_origin, aperture_par, manipulator_pos);
    */
  
    // Thanks to Keitaro Yamashita
    detx /= 50; dety /= 50;
    rotation *= M_PI / 180;
    
    int row = det_id - 1;

    // Nphotons = S [ADU] * G [e-/ADU] / (E [eV] * 3.65 [eV/e-]) according to the manual.
    // Thus, ADU/eV = 1/(3.65*G)

    //    ss << boost::format("q%d/adu_per_eV = %f\n") % det_id % (1.0 / (gain * 3.65));
    ss << boost::format("q%d/adu_per_eV = %f\n") % det_id % (1.0 / (0.1 * energy * 1000)); // Keitaro's 0.1 photon
    ss << boost::format("q%d/min_fs = %d\n") % det_id % 0;
    ss << boost::format("q%d/min_ss = %d\n") % det_id % (row * ysize);
    ss << boost::format("q%d/max_fs = %d\n") % det_id % 511;
    ss << boost::format("q%d/max_ss = %d\n") % det_id % ((row + 1) * ysize - 1);
    
    ss << boost::format("q%d/fs = %fx %+fy\n") % det_id % std::cos(rotation) % std::sin(rotation);
    ss << boost::format("q%d/ss = %fx %+fy\n") % det_id % -std::sin(rotation) % std::cos(rotation);
    ss << boost::format("q%d/corner_x = %f\n") % det_id % -detx;
    ss << boost::format("q%d/corner_y = %f\n\n") % det_id % dety;
  }
  return ss.str();
}

static void get_geom_h5(int runid, int taghi, int tag) {
  std::stringstream ss;
  float gain, detx, dety, detz, rotation;
//  float stage_dir, shift_weight, aperture_origin, aperture_par;
//  int manipulator_pos;
  double energy;
  
  printf("Panel geometry:\n");
  char det_name[256];
  for (int det_id = 1; det_id <= 8; det_id++) {
    snprintf(det_name, 256, det_name_template[det_temp_idx], det_id);
    
    ReadAbsGain(gain, det_name, bl, runid, taghi, tag);
    ReadDetPosX(detx, det_name, bl, runid, taghi, tag);
    ReadDetPosY(dety, det_name, bl, runid, taghi, tag);
    ReadDetPosZ(detz, det_name, bl, runid, taghi, tag);
    ReadDetRotationAngle(rotation, det_name, bl, runid, taghi, tag);
    
    ReadConfigOfPhotonEnergy(energy, bl, taghi, tag);
    printf(" panel %s gain %f pos (%f, %f, %f) rotation %f energy %f\n", det_name, gain, detx, dety, detz, rotation, energy);
  
    // Thanks to Keitaro Yamashita
    rotation *= M_PI / 180;
    
    int row = det_id - 1;
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
  printf("\n");

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

 static int add_image(double *buffer, int run, int taghi, int tag, double photon_energy) {
  int retno = 0;
  float buf_panel[xsize * ydatasize];
  float gain;
  
  char det_name[256];
  if (gains[0] == 0) { // Caching didn't improve performance
    //if (1) { 
    for (int det_id = 1; det_id <= 8; det_id++) {
      snprintf(det_name, 256, det_name_template[det_temp_idx], det_id);
      ReadAbsGain(gains[det_id], det_name, bl, run, taghi, tag);
    }
    gains[0] = 1;
  }

  /*  double photon_energy;
  ReadConfigOfPhotonEnergy(photon_energy, bl, taghi, tag); // returned in keV
  photon_energy *= 1000; // to eV */
  
  for (int det_id = 1; det_id <= 8; det_id++) {
    snprintf(det_name, 256, det_name_template[det_temp_idx], det_id);
    retno = ReadDetData(buf_panel, det_name, bl, run, taghi, tag);//, "calib");
    if (retno != 0) {
      printf("Tag %d not found.\n", tag);
      return 0;
    }

    gain = gains[det_id] / gains[1];
    gain *= gains[1] * 3.65 / 0.1 / photon_energy; // CHECKME: Keitaro's 0.1 photon unit
    //    printf("%f * %.10f ", buf_panel[0], gain);

    // Mark origin
    if (false) {
      for (int i = 0; i < 10; i++) {
	for (int j = 0; j < 10; j++) {
	  buf_panel[i * xsize + j] = det_id * 500;
	}
      }
    }
    
    int offset = (det_id - 1) * blocksize;
    for (int i = 0; i < blocksize; i++) {
      buffer[offset + i] += buf_panel[i] * gain;
    }
  }
  
  return 1;
}

