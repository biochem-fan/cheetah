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
#include "DataAccessUserAPI.h"
#include "hdf5.h"

const int NDET = 8;
const int parallel_size = 3; // MUST match dispatcher
const int bl = 3;
const int xsize = 512;
const int ysize = 1024; 
const int ydatasize = 1030;
const int blocksize = xsize * ysize;
const int buffersize = blocksize * NDET;
const int PD_ANY = -1, PD_LIGHT = 0, PD_DARK1 = 1, PD_DARK2 = 2;

char *det_name_template[30] = {"EventInfo_stor0%d", "MPCCD-8-2-001-%d", "EventInfo_stor1_0%d"};
char LLF_ID[50] = {};
char *LLF_ST4 = "BL3-ST4-MPCCD-octal";
char *LLF_ST3 = "BL3-ST3-MPCCD-octal";

// FIXME: make these local. Is Cheetah's main portion reentrant?
double buffer[buffersize] = {};
unsigned short averaged[buffersize] = {};
int det_temp_idx = -1;
float gains[9] = {};
char *streaders[NDET], *databufs[NDET];

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

// photon enery in eV
static bool get_image(double *buffer, int tag, double photon_energy) {
  int retno = 0;
  float buf_panel[xsize * ydatasize];
  float gain;

  for (int det_id = 0; det_id < NDET; det_id++) {
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
    gain *= gains[0] * 3.65 / 0.1 / photon_energy; // Keitaro's 0.1 photon unit

    // DEBUG: Mark origin
    if (false) {
      for (int i = 0; i < 10; i++) {
		  for (int j = 0; j < 10; j++) {
			  buf_panel[i * xsize + j] = det_id * 500;
		  }
      }
    }
    
    int offset = det_id * blocksize;
    for (int i = 0; i < blocksize; i++) {
      buffer[offset + i] = buf_panel[i] * gain;
    }
  }
  
  return true;
}

int main(int argc, char *argv[]) {
	printf("Cheetah for SACLA new offline API -- version 160105\n");
	printf(" by Takanori Nakane\n");
	printf(" This program is based on cheetah-sacla by Anton Barty.\n");
	int c, retno;
	
	// default values
	int runNumber = -1;
	char cheetahIni[4096] = {};
	int maxI_threshold = 10000;
	int startAt = 0;
	int stride = 2; // 30 Hz mode (60 / 30)
	double pd1_threshold = 0, pd2_threshold = 0;
	char *pd1_sensor_name = "xfel_bl_3_st_4_pd_laser_fitting_peak/voltage";
	char *pd2_sensor_name = "xfel_bl_3_st_4_pd_user_10_fitting_peak/voltage";
	int parallel_block = -1;
	int light_dark = PD_ANY;
	
	char outputH5[4096] = {};
	char station = 4;
	char *tagList_file = NULL;

	const struct option longopts[] = {
		{"ini", 1, NULL, 'i'},
		{"output", 1, NULL, 'o'},
		{"run", 1, NULL, 'r'},
		{"maxI", 1, NULL, 'm'},
		{"stride", 1, NULL, 11},
		{"station", 1, NULL, 12},
		{"pd1_thresh", 1, NULL, 13},
		{"pd2_thresh", 1, NULL, 14},
		{"type", 1, NULL, 15},
		{"list", 1, NULL, 16},
		{"pd1_name", 1, NULL, 17},
		{"pd2_name", 1, NULL, 18},
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
				if (parallel_block < -1 || parallel_block >= parallel_size) {
					printf("ERROR: wrong parallel_block.\n");
					return -1;
				}
			}
			break;
		case 16: // list
			if (tagList_file != NULL) {
				printf("ERROR: you cannot specify more than one tag list.\n");
				return -1;
			}
			tagList_file = strdup(optarg);
			printf("A tag list file was specified. maxI check was disabled.\n");
			maxI_threshold = -1;
			break;
		case 17: // pd1_name
			pd1_sensor_name = strdup(optarg); // small leak.
			break;
		case 18: // pd2_name
			pd2_sensor_name = strdup(optarg); // small leak.
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
	printf(" stride (--stride):            %d (default = 2; 30 Hz mode)\n", stride);
	printf(" station (--station):          %d (default = 4)\n", station);
	printf(" PD1 threshold (--pd1_thresh): %.3f (default = 0; ignore.)\n", pd1_threshold);
	printf(" PD2 threshold (--pd2_thresh): %.3f (default = 0; ignore.)\n", pd2_threshold);
	printf(" PD1 sensor name (--pd1_name): %s\n)", pd1_sensor_name);
	printf(" PD2 sensor name (--pd2_name): %s\n)", pd2_sensor_name);
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
    dims[0] = NDET * ysize;
    dims[1] = xsize;

	// get tag_hi and start
	int tag_hi, start, end;
	retno = sy_read_start_tagnumber(&tag_hi, &start, bl, runNumber);
	if (retno != 0) {
		printf("ERROR: Cannot read run %d.\n", runNumber);
		printf("If this run is before Nov 2014, please use the old API version.\n");
		return -1;
	}
	retno = sy_read_end_tagnumber(&tag_hi, &end, bl, runNumber);

	// How many dark frames?
	int numAll = (end - start) / stride + 1;
	printf("Run %d contains tag %d - %d (%d images), taghi = %d\n", runNumber, start, end, numAll, tag_hi);
	int *tagAll = (int*)malloc(sizeof(int) * numAll);
	for (int i = 0; i < numAll; i++) {
		tagAll[i] = start + i * stride;
	}

	std::vector<std::string> shutterAll;
	if (myReadSyncDataList(&shutterAll, "xfel_bl_3_shutter_1_open_valid/status", tag_hi, numAll, tagAll) != 0) {
		printf("Failed to get shutter status.\n");
		return -1;
	}

	int numDark = 0;
	for (int i = 0; i < numAll; i++) {
		if (atoi(shutterAll[i].c_str()) == 0) {
			numDark++;
		} else {
			break;
		}
	}
	printf("Number of dark frames: %d\n\n", numDark);
	startAt = (numDark + 1) * stride;
	free(tagAll);

	// find detector ID
	struct da_string_array *det_ids;
	int n_detid;

	printf("Detector configulation:\n");
	da_alloc_string_array(&det_ids);
	sy_read_detidlist(det_ids, bl, runNumber);
  
	da_getsize_string_array(&n_detid, det_ids);
	for (int i = 0; i < n_detid; i++) {
		char *detid;
		da_getstring_string_array(&detid, det_ids, i);
		printf(" detID #%d = %s\n", i, detid);
		if (strcmp(detid, "MPCCD-8-2-001-1")) {
			det_temp_idx = 1;
		} else if (strcmp(detid, "EventInfo_stor01")) {
			det_temp_idx = 0;
		}
		free(detid);    
	}
	if (det_temp_idx == -1) {
		printf("ERROR: Unknown detector ID.\n");
        cheetahExit(&cheetahGlobal);
        snprintf(message, 512, "Status=Error-BadDetID");
        cheetahGlobal.writeStatus(message);
		return -1;
	}
	da_destroy_string_array(&det_ids);
	
	start += startAt;
	int parallel_cnt = 0;
	std::vector<int> tagList;
	if (tagList_file == NULL) {
		int blockstart = start, blockend = end; // inclusive
		if (parallel_block != -1) { // block division
			int width = (end - start + 1) / parallel_size;
			blockstart = start + width * parallel_block;
			blockend = start + width * (parallel_block + 1) - 1;
			if (parallel_block == parallel_size - 1) { // last block
				blockend = end;
			}
		}
		printf("parallel: start %d end %d blockstart %d blockend %d\n", start, end, blockstart, blockend);
		for (int i = start; i <= end; i+= stride) {
			parallel_cnt++;
//			if (parallel_block == -1 || // no parallelization
//				parallel_cnt % parallel_size == parallel_block) {
			if (blockstart <= i && i <= blockend) {
				tagList.push_back(i);
			}
		}
	} else {
		FILE *fh = fopen(tagList_file, "r");
		free(tagList_file);

		if (fh == NULL) {
			printf("ERROR: Unable to open tagList.\n");
			return -1;
		}
		int i = 0;
		while (!feof(fh)) {
			fscanf(fh, "%d\n", &i);
			if (i < start || i > end) {
				printf("WARNING: tag %d does not belong to run %d. skipped.\n", i, runNumber);
				continue;
			}
 
			parallel_cnt++; // TODO: refactor and merge above
			if (parallel_block == -1 || // no parallelization
				parallel_cnt % parallel_size == parallel_block) {
				tagList.push_back(i);
			}
		}
		fclose(fh);
	}
//		tagList.clear(); tagList.push_back(121943650); // for debugging
	
	if (tagList.size() == 0) {
		printf("No images to process! Exiting...\n");
		cheetahExit(&cheetahGlobal);
		snprintf(message, 512, "Status=Error-NoImage");
		cheetahGlobal.writeStatus(message);
		return -1;
	}

	// for API 
	int *tagList_array = (int*)malloc(sizeof(int) * tagList.size());
	std::copy(tagList.begin(), tagList.end(), tagList_array);
	
	// Create storage readers and buffers
	// and get detector gains
	printf("Initializing reader and buffer\n");
	for (int det_id = 0; det_id < NDET; det_id++) {
		char det_name[256];
		snprintf(det_name, 256, det_name_template[det_temp_idx], det_id + 1);
		
		printf(" detector %s\n", det_name);
		retno = st_create_streader(&streaders[det_id], det_name, bl, 1, &runNumber);
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
		mp_read_absgain(&gains[det_id], databufs[det_id]);
	}

	// LLF values
	std::vector<std::string> maxIs;
	if (runNumber >= 355387 && runNumber <=355403) {
		// Broken runs in 2015-Jul Beamtime
		printf("WARNING: LLF ignored as the database is broken!\n");
	} else {
		struct da_string_array *llf;
		int n_llf;
		
		da_alloc_string_array(&llf);
		sy_read_statistics_detllf(llf, LLF_ID, bl, tag_hi, tagList.size(), tagList_array);
		
		da_getsize_string_array(&n_llf, llf);
		for (int i = 0; i < n_llf; i++) {
			char *val;
			da_getstring_string_array(&val, llf, i);
			maxIs.push_back(val);
			free(val);
		}
		da_destroy_string_array(&llf);
	}

	// Pulse energies (in keV)
	std::vector<std::string> pulse_energies;
    if (runNumber >=333661 && runNumber <= 333682) {
		// 19 May 2015: spectrometer broken! use config value instead
		printf("Using 7000 eV to fill in missing photon energies due to DB failure during run 333661-333682\n");
		for (unsigned int i = 0; i < tagList.size(); i++) {
			pulse_energies.push_back("7.0");
		}
	} else {
		retno = myReadSyncDataList(&pulse_energies, "xfel_bl_3_tc_spec_1/energy", tag_hi, tagList.size(), tagList_array);
		if (retno != 0) {
			printf("Failed to get photon_energy. Exiting...\n");
			cheetahExit(&cheetahGlobal);
			snprintf(message, 512, "Status=Error-PhotonEnergy");
			cheetahGlobal.writeStatus(message);
			return -1;
		}
	}

	std::vector<std::string> pd1_values, pd2_values, shutter;
	retno = myReadSyncDataList(&pd1_values, pd1_sensor_name, tag_hi, tagList.size(), tagList_array);
	if (retno != 0) {
		printf("WARNING: Failed to get %s.\n", pd1_sensor_name);
	}
	retno = myReadSyncDataList(&pd2_values, pd2_sensor_name, tag_hi, tagList.size(), tagList_array);
	if (retno != 0) {
		printf("WARNING: Failed to get %s.\n", pd2_sensor_name);
	}
	retno = myReadSyncDataList(&shutter, "xfel_bl_3_shutter_1_open_valid/status", tag_hi, tagList.size(), tagList_array);
	if (retno != 0) {
		printf("WARNING: Failed to get xfel_bl_3_shutter_1_open_valid/status.\n");
	}


	int processedTags = 0, LLFpassed = 0, tagSize = tagList.size(), frame_after_light = 0;
	for (int j = 0; j < tagSize; j++) {
		int tagID = tagList[j];
		int maxI = 0;
		if (runNumber >= 355387 && runNumber <=355403) {
			maxI = 100000; // Accept all
		} else {
			maxI = atoi(maxIs[j].c_str());
		}

		printf("tag %d shutter = %s\n", tagID, shutter[j].c_str());
		if (runNumber >= 358814 && runNumber <=358842) {
			// 2015 Oct: new run control GUI produces gaps in tag number
			if (atoi(shutter[j].c_str()) != 1) {
				printf("SHUTTER: tag %d rejected. shutter = %s\n", tagID, shutter[j].c_str());
				continue;
			}
		}

		double pd1_value = atof(pd1_values[j].c_str());
		double pd2_value = atof(pd2_values[j].c_str());
		double photon_energy; // in eV
		photon_energy = 1000 * atof(pulse_energies[j].c_str());

		bool light = true;
		if (pd1_threshold != 0 && 
			!(pd1_threshold > 0 && pd1_threshold <= pd1_value) &&
			!(pd1_threshold < 0 && -pd1_threshold > pd1_value)) light = false;
		if (pd2_threshold != 0 &&
			!(pd2_threshold > 0 && pd2_threshold <= pd2_value) &&
			!(pd2_threshold < 0 && -pd2_threshold > pd2_value)) light = false;
		if (light) frame_after_light = 0;
		else frame_after_light++;
		printf("Event %d: energy %f frame_after_light %d pd1_value %f pd2_value %f\n", tagID, photon_energy, frame_after_light, pd1_value, pd2_value);
		if ((light_dark == PD_DARK1 && frame_after_light != 1) ||
			(light_dark == PD_DARK2 && frame_after_light != 2) ||
			(light_dark == PD_LIGHT && frame_after_light != 0)) continue;

		processedTags++;

		printf("Event: %d (%d / %d (%.1f%%), Filter passed %d / %d (%.1f%%), Hits %ld (%.1f%%), maxI = %d, pd1_value = %.1f, pd2_value = %.3f\n",
			   tagID, (j + 1), tagSize, 100.0 * (j + 1) / tagSize, 
			   LLFpassed, processedTags, 100.0 * LLFpassed / processedTags,
			   cheetahGlobal.nhits, 100.0 * cheetahGlobal.nhits / processedTags,
			   maxI, pd1_value, pd2_value);
		if (maxI < maxI_threshold) {
			continue;
		}
		LLFpassed++;

		if (!get_image(buffer, tagID, photon_energy)) {
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

	free(tagList_array);
	
	time_t endT;
	time(&endT);
	double dif = difftime(endT,startT);
	std::cout << "time taken: " << dif << " seconds\n";
	std::cout << "Clean exit\n";
    return 0;
}

