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
int bl = 3;
const int xsize = 512;
const int ysize = 1024; 
const int ydatasize = 1030;
const int blocksize = xsize * ysize;
const int buffersize = blocksize * NDET;
const int PD_DARK_ANY = -2, PD_ANY = -1, PD_LIGHT = 0; // PD_DARKn = n afterwards

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
	printf("Cheetah for SACLA new offline API -- version 200123\n");
	printf(" by Takanori Nakane\n");
	printf(" This program is based on cheetah-sacla by Anton Barty.\n");
	int c, retno;
	
	// default values
	int runNumber = -1;
	char cheetahIni[4096] = {};
	double pd1_threshold = 0, pd2_threshold = 0, pd3_threshold = 0;
	char *pd1_sensor_name = "xfel_bl_3_st_4_pd_laser_fitting_peak/voltage";
	char *pd2_sensor_name = "xfel_bl_3_st_4_pd_user_10_fitting_peak/voltage";
	char *pd3_sensor_name = "xfel_bl_3_st_4_pd_user_10_fitting_peak/voltage"; // same as pd2 (dummy)
	int parallel_block = -1;
	int light_dark = PD_ANY;
	
	char outputH5[4096] = {};
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
		{"pd3_thresh", 1, NULL, 19},
		{"pd3_name", 1, NULL, 20},
                {"bl", 1, NULL, 21},
		{0, 0, NULL, 0}
	};

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
			printf("WARNING: LLF is no longer supported. Thus, maxI option is ignored.\n");
			break;
		case 11: //stride
			printf("WARNING: stride option is no longer supported; ignored.\n");
			break;
		case 12: // station
			printf("WARNING: LLF is no longer supported. Thus, station option is ignored.\n");
			break;
		case 13: // pd1_thresh
			pd1_threshold = atof(optarg);
			break;
		case 14: // pd2_thresh
			pd2_threshold = atof(optarg);
			break;
		case 19: // pd3_thresh
			pd3_threshold = atof(optarg);
			break;
		case 15: // type
			if (strcmp(optarg, "light") == 0) {
				light_dark = PD_LIGHT;
			} else if (strcmp(optarg, "dark") == 0) {
				light_dark = PD_DARK_ANY;
			} else if (strlen(optarg) == 5 &&
					   optarg[0] == 'd' && optarg[1] == 'a' &&
					   optarg[2] == 'r' && optarg[3] == 'k' ) { // darkN
				light_dark = optarg[4] - '0';
				if (light_dark < 1 || light_dark > 9) {
					printf("ERROR: wrong type.\n");
					return -1;
				}
			} else {
				parallel_block = atoi(optarg);
				if (parallel_block < -1 || parallel_block >= parallel_size) {
					printf("ERROR: wrong type or parallel_block.\n");
					return -1;
				}
				break;
			}
			break;
		case 16: // list
			if (tagList_file != NULL) {
				printf("ERROR: you cannot specify more than one tag list.\n");
				return -1;
			}
			tagList_file = strdup(optarg);
			printf("A tag list file was specified. maxI check was disabled.\n");
			break;
		case 17: // pd1_name
			pd1_sensor_name = strdup(optarg); // small leak.
			break;
		case 18: // pd2_name
			pd2_sensor_name = strdup(optarg); // small leak.
			break;
		case 20: // pd3_name
			pd3_sensor_name = strdup(optarg); // small leak.
			break;
                case 21: // bl
                        bl = atoi(optarg);
			if (bl != 2 && bl != 3) {
				printf("ERROR: beamline must be 2 or 3.\n");
				return -1;
			}
                        break;
		}
	}
	if (strnlen(outputH5, 4096) == 0) {
		snprintf(outputH5, 4096, "run%d.h5", runNumber);
	}

	printf("\nConfigurations:\n");
	printf(" runNumber (-r/--run):         %d\n", runNumber);
	printf(" cheetah ini file (-i/--ini):  %s\n", cheetahIni);
	printf(" output H5 file (-o/--output): %s (default = run######.h5)\n", outputH5);
        printf(" beamline (--bl):              %d (default = 3)\n", bl);
	printf(" PD1 threshold (--pd1_thresh): %.3f (default = 0; ignore.)\n", pd1_threshold);
	printf(" PD2 threshold (--pd2_thresh): %.3f (default = 0; ignore.)\n", pd2_threshold);
	printf(" PD3 threshold (--pd3_thresh): %.3f (default = 0; ignore.)\n", pd3_threshold);
	printf(" PD1 sensor name (--pd1_name): %s)\n", pd1_sensor_name);
	printf(" PD2 sensor name (--pd2_name): %s)\n", pd2_sensor_name);
	printf(" PD3 sensor name (--pd3_name): %s)\n", pd3_sensor_name);
	printf(" nFrame after light:           %d (default = -1; accept all image. -2; accept all dark images)\n", light_dark);
	printf(" parallel_block:               %d (default = -1; no parallelization)\n", parallel_block);

	if (runNumber < 0 || strlen(cheetahIni) == 0) {
		printf("Wrong argument! \nUsage: cheetah-sacla-api -i cheetah.ini -r runNumber\n");
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
	printf("\n");

	hsize_t dims[2];
	dims[0] = NDET * ysize;
	dims[1] = xsize;

	// get tag_hi and start
	int tag_hi, tag_start, tag_end;
	retno = sy_read_start_tagnumber(&tag_hi, &tag_start, bl, runNumber);
	if (retno != 0) {
		printf("ERROR: Cannot read run %d.\n", runNumber);
		printf("If this run is before Nov 2014, please use the old API version.\n");
		return -1;
	}
	retno = sy_read_end_tagnumber(&tag_hi, &tag_end, bl, runNumber);

	// How many dark frames?
	struct da_int_array *tagbuf;
	int numAll;
	da_alloc_int_array(&tagbuf, 0, NULL);
	sy_read_taglist_byrun(tagbuf, bl, runNumber);
	da_getsize_int_array(&numAll, tagbuf);
	printf("Run %d contains tag %d (inclusive) to %d (exclusive), thus %d images\n", runNumber, tag_start, tag_end, numAll);
	int *tagAll = (int*)malloc(sizeof(int) * numAll);
	for (int i = 0; i < numAll; i++) {
		da_getint_int_array(tagAll + i, tagbuf, i);
	}
	da_destroy_int_array(&tagbuf);

	// See prepare-cheetah-sacla-api2.py for details.
	bool workaround_18feb = (bl == 2) && (runNumber >= 32348) && (runNumber < 81550);
        bool use_close_status = (bl == 2) && (runNumber >= 81550);
	std::vector<std::string> shutterOpen, shutterClose;
	if (workaround_18feb) {
		if (myReadSyncDataList(&shutterOpen, "xfel_bl_2_st_3_bm_1_pd/charge", tag_hi, numAll, tagAll) != 0) {
			printf("Failed to get xfel_bl_2_st_3_bm_1_pd/charge.\n");
			return -1;
		}
	} else if (bl == 3) {
		if (myReadSyncDataList(&shutterOpen, "xfel_bl_3_shutter_1_open_valid/status", tag_hi, numAll, tagAll) != 0) {
			printf("Failed to get xfel_bl_3_shutter_1_open_valid/status.\n");
			return -1;
		}
        	if (use_close_status && myReadSyncDataList(&shutterClose, "xfel_bl_3_shutter_1_close_valid/status", tag_hi, numAll, tagAll) != 0) {
			printf("Failed to get xfel_bl_3_shutter_1_close_valid/status.\n");
			return -1;
		}
	} else if (bl == 2) {
		if (myReadSyncDataList(&shutterOpen, "xfel_bl_2_shutter_1_open_valid/status", tag_hi, numAll, tagAll) != 0) {
			printf("Failed to get xfel_bl_2_shutter_1_open_valid/status.\n");
			return -1;
		}
		if (use_close_status && myReadSyncDataList(&shutterClose, "xfel_bl_2_shutter_1_close_valid/status", tag_hi, numAll, tagAll) != 0) {
			printf("Failed to get xfel_bl_2_shutter_1_close_valid/status.\n");
			return -1;
		}
	}

	if (workaround_18feb) printf("Warning: applyied workaround for unreliable shutter status since 18 feb.\n");
	if (use_close_status) printf("Both shutter open and shutter close fields are used to find non-exposed images.\n");
	int numDark = 0;
	for (int i = 0; i < numAll; i++) {
		tag_start = tagAll[i];
		// xfel_bl_2_st_3_bm_1_pd/charge might contain "saturated" or "not-converged".
		if ((shutterOpen[i] != "saturated" && (shutterOpen[i] == "not-converged" || atoi(shutterOpen[i].c_str()) == 0))
		    && (!use_close_status || atoi(shutterClose[i].c_str()) == 1)) {
			numDark++;
		} else {
			break;
		}
	}
	printf("Number of dark frames: %d\nThus, exposed images start from %d\n\n", numDark, tag_start);

	int parallel_cnt = 0;
	std::vector<int> tagList;
	if (tagList_file == NULL) {
		int blockstart = numDark, blockend = numAll - 1; // inclusive
		if (parallel_block != -1) { // block division
			int width = (numAll - numDark + 1) / parallel_size;
			blockstart = numDark + width * parallel_block;
			blockend = numDark + width * (parallel_block + 1) - 1;
			if (parallel_block == parallel_size - 1) { // last block
				blockend = numAll - 1;
			}
		}
		printf("parallel: start %d end %d blockstart %d blockend %d\n", tagAll[0], tagAll[numAll - 1], tagAll[blockstart], tagAll[blockend]);
		for (int i = blockstart; i <= blockend; i++) {
			tagList.push_back(tagAll[i]);
		
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
			if (i < tag_start || i >= tag_end) {
				printf("WARNING: tag %d does not belong to run %d. skipped.\n", i, runNumber);
				continue;
			}
 
			parallel_cnt++; // TODO: refactor and merge above, use block parallelization
			if (parallel_block == -1 || // no parallelization
				parallel_cnt % parallel_size == parallel_block) {
				tagList.push_back(i);
			}
		}
		fclose(fh);
	}
        free(tagAll);
	
	if (tagList.size() == 0) {
		printf("No images to process! Exiting...\n");
		cheetahExit(&cheetahGlobal);
		snprintf(message, 512, "Status=Error-NoImage");
		cheetahGlobal.writeStatus(message);
		return -1;
	}
	printf("\n");

	// for API 
	int *tagList_array = (int*)malloc(sizeof(int) * tagList.size());
	std::copy(tagList.begin(), tagList.end(), tagList_array);
	
	// find detector ID
	struct da_string_array *det_ids;
	int n_detid;
	char det_template[256] = {};

	printf("Detector configulation:\n");
	da_alloc_string_array(&det_ids);
	sy_read_detidlist(det_ids, bl, runNumber);

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
	da_destroy_string_array(&det_ids);

	if (det_template[0] == 0) {
		printf("ERROR: Unknown or non-supported detector ID.\n");
	        cheetahExit(&cheetahGlobal);
	        snprintf(message, 512, "Status=Error-BadDetID");
		cheetahGlobal.writeStatus(message);
		return -1;
	}
	printf("\n");
	
	// Create storage readers and buffers
	// and get detector gains
	printf("Initializing reader and buffer\n");
	for (int det_id = 0; det_id < NDET; det_id++) {
		char det_name[256];
		strncpy(det_name, det_template, 255);
		det_name[strlen(det_name) - 1] = '0' + det_id + 1;
		
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
		uint32_t tagid = tag_start;
		retno = st_collect_data(databufs[det_id], streaders[det_id], &tagid);
		if (retno != 0) {
			printf("Failed to collect data for %s.\n", det_name);
			return -1;
		}
		mp_read_absgain(&gains[det_id], databufs[det_id]);
	}
	printf("\n");

	// Pulse energies (in keV)
	double config_photon_energy;
	sy_read_config_photonenergy(&config_photon_energy, bl, runNumber);
	
	std::vector<std::string> pulse_energies;
	if (runNumber >=333661 && runNumber <= 333682) {
		// 19 May 2015: spectrometer broken! use config value instead
		printf("Using 7000 eV to fill in missing photon energies due to DB failure during run 333661-333682\n");
		for (unsigned int i = 0; i < tagList.size(); i++) {
			pulse_energies.push_back("7.0");
		}
	} else {
                if (bl == 3) {
			retno = myReadSyncDataList(&pulse_energies, "xfel_bl_3_tc_spec_1/energy", tag_hi, tagList.size(), tagList_array);
		} else if (bl == 2) {
			retno = myReadSyncDataList(&pulse_energies, "xfel_bl_2_tc_spec_1/energy", tag_hi, tagList.size(), tagList_array);
		}
		if (retno != 0) {
			printf("Failed to get photon_energy. Exiting...\n");
			cheetahExit(&cheetahGlobal);
			snprintf(message, 512, "Status=Error-PhotonEnergy");
			cheetahGlobal.writeStatus(message);
			return -1;
		}
	}

	std::vector<std::string> pd1_values, pd2_values, pd3_values, shutter;
	retno = myReadSyncDataList(&pd1_values, pd1_sensor_name, tag_hi, tagList.size(), tagList_array);
	if (retno != 0) {
		printf("WARNING: Failed to get %s.\n", pd1_sensor_name);
	}
	retno = myReadSyncDataList(&pd2_values, pd2_sensor_name, tag_hi, tagList.size(), tagList_array);
	if (retno != 0) {
		printf("WARNING: Failed to get %s.\n", pd2_sensor_name);
	}
	retno = myReadSyncDataList(&pd3_values, pd3_sensor_name, tag_hi, tagList.size(), tagList_array);
	if (retno != 0) {
		printf("WARNING: Failed to get %s.\n", pd3_sensor_name);
	}
	if (bl == 3) {
		retno = myReadSyncDataList(&shutter, "xfel_bl_3_shutter_1_open_valid/status", tag_hi, tagList.size(), tagList_array);
	} else if (bl == 2) {
		retno = myReadSyncDataList(&shutter, "xfel_bl_2_shutter_1_open_valid/status", tag_hi, tagList.size(), tagList_array);
	}
	if (retno != 0) {
		printf("WARNING: Failed to get shutter status.\n");
	}

	int processedTags = 0, LLFpassed = 0, tagSize = tagList.size(), frame_after_light = 9999;
	for (int j = 0; j < tagSize; j++) {
		int tagID = tagList[j];

		printf("tag %d shutter = %s\n", tagID, shutter[j].c_str());
		if (runNumber >= 358814 && runNumber <=358842) {
			// 2015 Oct: new run control GUI produces gaps in tag number
			//           attempts to read such tags cause crash, so we have to skip.
			// For other runs, reading shutter-closed images does not cause any harm.
			// Actually, some runs in 2018 Feb have unreliable shuter status, so
			// we should NOT skip images based on shutter status!

			if (atoi(shutter[j].c_str()) != 1) {
				printf("SHUTTER: tag %d rejected. shutter = %s\n", tagID, shutter[j].c_str());
				continue;
			}
		}

		double pd1_value = atof(pd1_values[j].c_str());
		double pd2_value = atof(pd2_values[j].c_str());
		double pd3_value = atof(pd3_values[j].c_str());
		double photon_energy; // in eV
		photon_energy = 1000 * atof(pulse_energies[j].c_str());
		if (photon_energy == 0) {
			printf("WARNING: The wavelength from the inline spectrometer for tag %d is not available\n", tagID);
			if (config_photon_energy < 5000 || config_photon_energy > 14000) {
				printf("         The accelerator config value also looks broken; assumed 7 keV as a last resort.\n");
				printf("         The scale factor for this image can be wrong!!\n");
				photon_energy = 7000;
			} else {
				printf("         Used the accelerator config value (%f eV) instead.\n", config_photon_energy);
				photon_energy = config_photon_energy;
			}
		}

		bool light = true;
		if (pd1_threshold != 0 && 
			!(pd1_threshold > 0 && pd1_threshold <= pd1_value) &&
			!(pd1_threshold < 0 && -pd1_threshold > pd1_value)) light = false;
		if (pd2_threshold != 0 &&
			!(pd2_threshold > 0 && pd2_threshold <= pd2_value) &&
			!(pd2_threshold < 0 && -pd2_threshold > pd2_value)) light = false;
		if (pd3_threshold != 0 &&
			!(pd3_threshold > 0 && pd3_threshold <= pd3_value) &&
			!(pd3_threshold < 0 && -pd3_threshold > pd3_value)) light = false;
		if (light) {
			frame_after_light = 0;
		} else {
			frame_after_light++;
		}
//		printf("Event %d: energy %f frame_after_light %d pd1_value %f pd2_value %f pd3_value %f\n", tagID, photon_energy, frame_after_light, pd1_value, pd2_value, pd3_value);
		if ((light_dark >= 0 && frame_after_light != light_dark) ||
			(light_dark == PD_DARK_ANY && frame_after_light == 0)) continue;

		processedTags++;

		printf("Event: %d (%d / %d (%.1f%%), Filter passed %d / %d (%.1f%%), Hits %ld (%.1f%%), pd1_value = %.1f, pd2_value = %.3f, pd3_value = %.3f\n",
			   tagID, (j + 1), tagSize, 100.0 * (j + 1) / tagSize, 
			   LLFpassed, processedTags, 100.0 * LLFpassed / processedTags,
			   cheetahGlobal.nhits, 100.0 * cheetahGlobal.nhits / processedTags,
			   pd1_value, pd2_value, pd3_value);
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

