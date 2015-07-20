//
//  main.cpp
//  cheetah-sacla-api
//
//  By Takanori Nakane, based on cheetah-sacla-hdf5
//
//  Created by Anton Barty on 20/1/14.
//  Copyright (c) 2014 Anton Barty. All rights reserved.
//

#include <iostream>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include "cheetah.h"
#include "SaclaDataAccessUserAPI.h"
#include "hdf5.h"

const int parallel_size = 3; // MUST match dispatcher
const int bl = 3;
const int xsize = 512;
const int ysize = 1024; 
const int ydatasize = 1030;
const int blocksize = xsize * ysize;
const int buffersize = blocksize * 8;
const int PD_ANY = -1, PD_LIGHT = 0, PD_DARK1 = 1, PD_DARK2 = 2;

char *det_name_template[30] = {"EventInfo_stor0%d", "MPCCD-8-2-001-%d", "EventInfo_stor1_0%d"};
char LLF_ID[50] = {};
char *LLF_ST4 = "BL3-ST4-MPCCD-octal";
char *LLF_ST3 = "BL3-ST3-MPCCD-octal";

// FIXME: make these local. Is Cheetah's main portion reentrant?
double buffer[buffersize] = {};
unsigned short averaged[buffersize] = {};
float posx[buffersize] = {}, posy[buffersize] = {}, posz[buffersize] = {};
int det_temp_idx = -1;
float gains[9] = {};

// photon enery in eV
static bool get_image(double *buffer, int run, int taghi, int tag, double photon_energy) {
  int retno = 0;
  float buf_panel[xsize * ydatasize];
  float gain, gain_panel1;
  char det_name[256];

  if (gains[0] == 0) { // Caching didn't improve performance
    for (int det_id = 1; det_id <= 8; det_id++) {
		snprintf(det_name, 256, det_name_template[det_temp_idx], det_id);
		ReadAbsGain(gains[det_id], det_name, bl, run, taghi, tag);
		printf("Detector absolute gain for panel %s = %f\n", det_name, gains[det_id]);
    }
    gains[0] = 1;
  }

  for (int det_id = 1; det_id <= 8; det_id++) {
    snprintf(det_name, 256, det_name_template[det_temp_idx], det_id);
    retno = ReadDetData(buf_panel, det_name, bl, run, taghi, tag, "calib");
    if (retno != 0) {
      printf("Tag %d not found.\n", tag);
      return 0;
    }

    gain = gains[det_id] / gains[1];
    gain *= gains[1] * 3.65 / 0.1 / photon_energy; // Keitaro's 0.1 photon unit

    // DEBUG: Mark origin
    if (false) {
      for (int i = 0; i < 10; i++) {
		  for (int j = 0; j < 10; j++) {
			  buf_panel[i * xsize + j] = det_id * 500;
		  }
      }
    }
    
    int offset = (det_id - 1) * blocksize;
    for (int i = 0; i < blocksize; i++) {
      buffer[offset + i] = buf_panel[i] * gain;
    }
  }
  
  return true;
}

int main(int argc, char *argv[]) {
	printf("Cheetah using SACLA API by Takanori Nakane\n");
	printf("This program is based on cheetah-sacla by Anton Barty\n");

	int c, retno;
	
	// default values
	int runNumber = -1;
	char cheetahIni[4096] = {};
	int maxI_threshold = 10000;
	int startAt = 300;
	int stride = 2; // 30Hz mode (60 / 30)
	double pd1_threshold = 0, pd2_threshold = 0;
	int parallel_block = -1;
	int light_dark = PD_ANY;
	
	char outputH5[4096] = {};
	char station = 4;

	const struct option longopts[] = {
		{"ini", 1, NULL, 'i'},
		{"output", 1, NULL, 'o'},
		{"run", 1, NULL, 'r'},
		{"maxI", 1, NULL, 'm'},
		{"start_at", 1, NULL, 10},
		{"stride", 1, NULL, 11},
		{"station", 1, NULL, 12},
		{"pd1_thresh", 1, NULL, 13},
		{"pd2_thresh", 1, NULL, 14},
		{"type", 1, NULL, 15},
		{0, 0, NULL, 0}
	};

	int tmp;
	while ((c = getopt_long(argc, argv, "i:r:m:o:", longopts, NULL)) != -1) {
		switch(c) {
		case 'i':
			strncpy(cheetahIni, optarg, 4096);
			break;
		case 'o':
			strncpy(outputH5, optarg, 4096);
			break;
		case 'r':
			runNumber = atoi(optarg);
			break;
		case 'm':
			maxI_threshold = atoi(optarg);
			break;
		case 10: // start-at
			startAt = atoi(optarg);
			break;
		case 11: // stride
			tmp = atoi(optarg);
			if (tmp <= 0) {
				printf("ERROR: Stride must be a positive integer. \n");
				return -1;
			}
			stride = tmp;
			break;
		case 12: // station
			station = atoi(optarg);
			break;
		case 13: // pd1_thresh
			pd1_threshold = atof(optarg);
			break;
		case 14: // pd2_thresh
			pd2_threshold = atof(optarg);
			break;
		case 15: // type
			if (strcmp(optarg, "light") == 0) {
				light_dark = PD_LIGHT;
			} else if (strcmp(optarg, "dark1") == 0) {
				light_dark = PD_DARK1;
			} else if (strcmp(optarg, "dark2") == 0) {
				light_dark = PD_DARK2;
			} else {
				parallel_block = atoi(optarg);
				if (parallel_block < 0 || parallel_block >= parallel_size) {
					printf("ERROR: wrong parallel_block.\n");
					return -1;
				}
			}
			break;
		}
	}
	if (strnlen(outputH5, 4096) == 0) {
		snprintf(outputH5, 4096, "run%d.h5", runNumber);
	}
	if (station == 3) {
		strncpy(LLF_ID, LLF_ST3, 50);
	} else if (station == 4) {
		strncpy(LLF_ID, LLF_ST4, 50);
	} else {
		printf("ERROR! Station must be 3 or 4\n");
		return -1;
	}

	printf("\nConfigurations:\n");
	printf(" runNumber (-r/--run):         %d\n", runNumber);
	printf(" cheetah ini file (-i/--ini):  %s\n", cheetahIni);
	printf(" output H5 file (-o/--output): %s (default = run######.h5)\n", outputH5);
	printf(" maxI threshold (-m/--maxI):   %d (default = 10000)\n", maxI_threshold);
	printf(" start frame (--start_at):     %d (default = 300; to skip 150 dark frames)\n", startAt);
	printf(" stride (--stride):            %d (default = 2; 30Hz mode)\n", stride);
	printf(" station (--station):          %d (default = 4)\n", station);
	printf(" PD1 threshold (--pd1_thresh): %.3f (default = 0; ignore.)\n", pd1_threshold);
	printf(" PD2 threshold (--pd2_thresh): %.3f (default = 0; ignore.)\n", pd2_threshold);
	printf(" nFrame after light:           %d (default = -1; any)\n", light_dark);
	printf(" parallel_block:               %d (default = -1; no parallelization)\n", parallel_block);


	if (runNumber < 0 || strlen(cheetahIni) == 0) {
		printf("Wrong argument! \nUsage: cheetah-sacla-api -i cheetah.ini -r runNumber [-m maxI]\n");
		return -1;
	}
	
	/*
	 *	Initialise Cheetah
	 */
	printf("\nSetting up Cheetah...\n");
	static uint32_t ntriggers = 0;
	static long frameNumber = 0;
	static cGlobal cheetahGlobal;
	char message[512];
	static time_t startT = 0;
	time(&startT);
    strcpy(cheetahGlobal.configFile, cheetahIni);
	cheetahInit(&cheetahGlobal);
	cheetahGlobal.runNumber = runNumber;
	strncpy(cheetahGlobal.cxiFilename, outputH5, MAX_FILENAME_LENGTH);

    hsize_t dims[2];
    dims[0] = 8 * ysize;
    dims[1] = xsize;
    
	int tag_hi, start, end;
	ReadStartTagNumber(tag_hi, start, bl, runNumber);
	retno = ReadEndTagNumber(tag_hi, end, bl, runNumber);
	printf("tag_hi = %d start = %d end = %d retno = %d\n", tag_hi, start, end, retno);
	
	std::vector<std::string> det_ids;
	ReadDetIDList(&det_ids, bl, runNumber);  
	
	for (unsigned int i = 0; i < det_ids.size(); i++) {
		if (det_ids[i] == "MPCCD-8-2-001-1") {
			det_temp_idx = 1;
			break;
		} else if (det_ids[i] == "EventInfo_stor01") {
			det_temp_idx = 0;
			break;
		} /*else if (det_ids[i] == "EventInfo_stor1_01") {
      det_temp_idx = 2;
      break;
      }*/
	}
	if (det_temp_idx == -1) {
		printf("ERROR: Unknown detector ID %s.\n", det_ids[0].c_str());
        cheetahExit(&cheetahGlobal);
        snprintf(message, 512, "Status=Error-BadDetID");
        cheetahGlobal.writeStatus(message);
		return -1;
	}
	
	start += startAt;
	std::vector<int> tagList;
	int parallel_cnt = 0;
	for (int i = start; i <= end; i+= stride) {
		parallel_cnt++;
		if (parallel_block == -1 || // no parallelization
			parallel_cnt % parallel_size == parallel_block)
			tagList.push_back(i);
	}
//		tagList.clear(); tagList.push_back(121943650); // for debugging
	
	std::vector<std::string> maxIs;
	if (tagList.size() == 0) {
		printf("No images to process! Exiting...\n");
		cheetahExit(&cheetahGlobal);
		snprintf(message, 512, "Status=Error-NoImage");
		cheetahGlobal.writeStatus(message);
		return -1;
	}
	ReadStatisticsOfDetLLF(&maxIs, LLF_ID, bl, tag_hi, tagList);
	
	std::vector<std::string> pulse_energies;
    if (runNumber >=333661 && runNumber <= 333682) {
		// 19 May 2015: spectrometer broken! use config value instead
		printf("Using 7000 eV to fill in missing photon energies due to DB failure during run 333661-333682\n");
		for (int i = 0; i < tagList.size(); i++) {
			pulse_energies.push_back("7.0");
		}
	} else {
		retno = ReadSyncDataList(&pulse_energies, "xfel_bl_3_tc_spec_1/energy", tag_hi, tagList);
		if (retno != 0) {
			printf("Failed to get photon_energy. Exiting...\n");
			cheetahExit(&cheetahGlobal);
			snprintf(message, 512, "Status=Error-PhotonEnergy");
			cheetahGlobal.writeStatus(message);
			return -1;
		}
	}

	std::vector<std::string> pd_laser, pd_user2;
	retno = ReadSyncDataList(&pd_laser, "xfel_bl_3_st_4_pd_laser_fitting_peak/voltage", tag_hi, tagList);
	if (retno != 0) {
		printf("WARNING: Failed to get xfel_bl_3_st_4_pd_laser_fitting_peak/voltage.\n");
	}
	retno = ReadSyncDataList(&pd_user2, "xfel_bl_3_st_4_pd_user_10_fitting_peak/voltage", tag_hi, tagList);
	if (retno != 0) {
		printf("WARNING: Failed to get xfel_bl_3_st_4_pd_user_10_fitting_peak/voltage.\n");
	}

	int processedTags = 0, LLFpassed = 0, tagSize = tagList.size(), frame_after_light = 0;
	for (int j = 0; j < tagSize; j++) {
		int tagID = tagList[j];
		int maxI = atoi(maxIs[j].c_str());
		double pd_laser_val = atof(pd_laser[j].c_str());
		double pd_user2_val = atof(pd_user2[j].c_str());
		double photon_energy; // in eV
		photon_energy = 1000 * atof(pulse_energies[j].c_str());

		bool light = true;
		if (pd1_threshold != 0 && 
			!(pd1_threshold > 0 && pd1_threshold <= pd_laser_val) &&
			!(pd1_threshold < 0 && -pd1_threshold > pd_laser_val)) light = false;
		if (pd2_threshold != 0 &&
			!(pd2_threshold > 0 && pd2_threshold <= pd_user2_val) &&
			!(pd2_threshold < 0 && -pd2_threshold > pd_user2_val)) light = false;
		if (light) frame_after_light = 0;
		else frame_after_light++;
		printf("Event %d: energy %f frame_after_light %d pd_user2_val %f\n", tagID, photon_energy, frame_after_light, pd_user2_val);
		if ((light_dark == PD_DARK1 && frame_after_light != 1) ||
			(light_dark == PD_DARK2 && frame_after_light != 2) ||
			(light_dark == PD_LIGHT && frame_after_light != 0)) continue;

		processedTags++;

		printf("Event: %d (%d / %d (%.1f%%), LLF passed %d / %d (%.1f%%), Hits %ld (%.1f%%), maxI = %d, PD_laser = %.1f, PD_user = %.3f\n",
			   tagID, (j + 1), tagSize, 100.0 * (j + 1) / tagSize, 
			   LLFpassed, processedTags, 100.0 * LLFpassed / processedTags,
			   cheetahGlobal.nhits, 100.0 * cheetahGlobal.nhits / processedTags,
			   maxI, pd_laser_val, pd_user2_val);
		if (maxI < maxI_threshold) {
			continue;
		}
		LLFpassed++;

		if (!get_image(buffer, runNumber, tag_hi, tagID, photon_energy)) {
			continue; // image not available
		}
		frameNumber++;
		
		/*
		 *  Cheetah: Calculate time beteeen processing of data frames
		 */
		time_t	tnow;
		double	dtime, datarate;
		time(&tnow);		
		dtime = difftime(tnow, cheetahGlobal.tlast);
		if(dtime > 1.) {
			datarate = (frameNumber - cheetahGlobal.lastTimingFrame)/dtime;
			cheetahGlobal.lastTimingFrame = frameNumber;
			time(&cheetahGlobal.tlast);
			cheetahGlobal.datarate = datarate;
		}

		snprintf(message, 512, "Total=%d,Processed=%d,LLFpassed=%d,Hits=%ld,Status=Hitfinding",
				 tagSize, (j + 1), LLFpassed, cheetahGlobal.nhits); 
		if (processedTags % 5 == 0) {
			cheetahGlobal.writeStatus(message);
		}
            			
		/*
		 *	Cheetah: Create a new eventData structure in which to place all information
		 */
		cEventData	*eventData;
		eventData = cheetahNewEvent(&cheetahGlobal);
		ntriggers++;
		
		eventData->pulnixFail = 1;
		eventData->specFail = 1;
		eventData->frameNumber = tagID;
		eventData->runNumber = runNumber;
		eventData->nPeaks = 0;
		eventData->pumpLaserCode = 0;
		eventData->pumpLaserDelay = 0;
		eventData->photonEnergyeV = photon_energy; // in eV
		eventData->wavelengthA = 12398 / eventData->photonEnergyeV; // 4.1357E-15 * 2.9979E8 * 1E10 / eV (A)
		eventData->pGlobal = &cheetahGlobal;
		eventData->fiducial = tagID; // must be unique
		
		int detID = 0;
		long pix_nn = cheetahGlobal.detector[detID].pix_nn;

//        int underflow = 0, overflow = 0;
		for(long ii = 0; ii < pix_nn; ii++) {
			long tmp = lrint(buffer[ii]);
			if (tmp < 0) {
//				underflow++; printf("%ld ", tmp);
				tmp = 0; 
			} else if (tmp > USHRT_MAX) {
//				overflow++;
				tmp = USHRT_MAX; 
			}
			eventData->detector[detID].data_raw16[ii] = (uint16_t)tmp;
		}
//		printf("#underflow = %d, #overflow = %d\n", underflow, overflow);
		
		cheetahProcessEventMultithreaded(&cheetahGlobal, eventData);
	}

	cheetahExit(&cheetahGlobal);
	snprintf(message, 512, "Total=%d,Processed=%d,LLFpassed=%d,Hits=%ld,Status=Finished",
			 tagSize, tagSize, LLFpassed, cheetahGlobal.nhits); 
	cheetahGlobal.writeStatus(message); // Overwrite "Status: Finished"
	
	time_t endT;
	time(&endT);
	double dif = difftime(endT,startT);
	std::cout << "time taken: " << dif << " seconds\n";
	std::cout << "Clean exit\n";
    return 0;
}

