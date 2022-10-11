//
//  cheetah-sacla-online-api
//
//  Originally created by Anton Barty on 20/1/14.
//  Copyright (c) 2014 Anton Barty. All rights reserved.
//
//  Modified by Takanori Nakane to read data from Online API
//
//  gaincal file was not used for the sake of CPU cache efficiency.

// TODO: apply gain correction (necessary??)

#include <iostream>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#include <pthread.h>
#include <time.h>
#include <map>

#include "cheetah.h"
#include "OnlineUserAPI.h"

void ol_initialize_dummy(char* filename);

#define XSIZE 512
#define YSIZE 1024
#define PANELSIZE (XSIZE * YSIZE)
#define IMGSIZE (PANELSIZE * 8)
#define SINGLE_THREAD false

#define NDET 10 // Normally 8
#define PRIMARY_DET 0 // Normally 0. 4 for xu-bl2-st3-opcon2
int TAG_INCREMENT = 2; // 2 for 30 Hz, 1 for API stub or 60 Hz

const struct timespec SHORT_WAIT = {0, 1000000}; // 1E6 = msec
const struct timespec LONG_WAIT = {0, 5000000}; // 25E6 for testing

typedef struct {
	float buf[IMGSIZE];
	int run;
	int error;
} Image_Data;

typedef struct {
	int detID;
	int socketID;
	pthread_mutex_t *mutex_map;
	std::map <int, Image_Data*> *images;
	int *cur_tags; // tag to look at NEXT
	float photon_energy;
} Thread_Data;

void* thread(void *thread_data) {
	// Initialize thread local variables

	int det_id = ((Thread_Data*)thread_data)->detID;
	int socket_id = ((Thread_Data*)thread_data)->socketID;
	pthread_mutex_t *mutex_map = ((Thread_Data*)thread_data)->mutex_map;
	std::map <int, Image_Data*> *images = ((Thread_Data*)thread_data)->images;
	int *cur_tags = ((Thread_Data*)thread_data)->cur_tags;
	float photon_energy = ((Thread_Data*)thread_data)->photon_energy;
	free(thread_data);

	printf("Child%d: created. socket = %d\n", det_id, socket_id);
	
	int datasize, worksize;
	ol_get_data_size(&datasize, &worksize, socket_id);
	printf("API: ol_get_data_size for %d returned dataSize = %d, workSize = %d\n", det_id, datasize, worksize);
	
	char *pDataStBuf = (char*)malloc(datasize);
	char *pWorkBuf = (char*)malloc(worksize);

	// Get the tag to start at
	
	if (det_id == PRIMARY_DET) {
		int actual_tag = -1;
		int err = ol_collect_det_data(pDataStBuf, pWorkBuf, &actual_tag, socket_id, -1, datasize, worksize);
		if (err < 0) {
			printf("ERROR: Failed to read the first image.\n");
			exit(-1);
		}
		printf("Start tag = %d\n", actual_tag);

		for (int i = 0; i < NDET; i++) {
			cur_tags[i] = actual_tag;
		}
	} else {
		while (cur_tags[det_id] == -1) {
			nanosleep(&SHORT_WAIT, NULL);
		}
	}

	// Main loop
	while (true) {
		std::map<int, Image_Data*>::const_iterator it;
		Image_Data *img = NULL;
		int wanted_tag = cur_tags[det_id];

		for (int ntry = 0; ntry < 200; ntry++) { // TODO: Check if this is enough
			pthread_mutex_lock(mutex_map);
			it = images->find(wanted_tag);
			if (it == images->end()) {
				if (det_id == PRIMARY_DET) {
					// - Only thread 0 allocates new buffer
					// - New buffer is allocated regardless of the tag status
					// This ensures that ALL even tags are enumerated.

					img = (Image_Data*)calloc(1, sizeof(Image_Data));
					img->error = false;
					(*images)[wanted_tag] = img;
					printf("Child%d: Allocated image %d\n", det_id, wanted_tag);
				}
			} else {
				img = it->second;
			}
			pthread_mutex_unlock(mutex_map);
			
			if (img == NULL) {
				// wait till the buffer is allocated
				nanosleep(&SHORT_WAIT, NULL);
                // printf("Child%d: Waiting for %d\n", det_id, wanted_tag);
			} else {
				break;
			}
		}
		
		if (img == NULL) {
			printf("Child%d: tag %d never get allocated. skipped.\n", det_id, wanted_tag);
		} else {
			int actual_tag = -1;
			int err = ol_collect_det_data(pDataStBuf, pWorkBuf, &actual_tag, socket_id, wanted_tag, datasize, worksize);
			if (err < 0 || actual_tag != wanted_tag) {
				printf("Child%d: could not get tag %d. skipped. error code %d\n", det_id, wanted_tag, err);
				img->error = true;
				exit(-1);
				// 151005: Exit right away; probably MPCCD server has been restarted.
				// Of course this might happen on our side, for example, if we come too late (TAGDATAGONE -10000).
				// But still it is worth restarting. 
			} else {
//				printf("Child%d: Got image %d\n", det_id, wanted_tag);
				if (det_id == PRIMARY_DET) {
					ol_read_run_num(&img->run, pDataStBuf, 0);
				}

				float *detData;
				ol_read_det_data(&detData, pDataStBuf, 0);
				float gain;
				ol_read_abs_gain(&gain, pDataStBuf, 0);
				gain *= 3.65 / 0.1 / photon_energy;
//				gain = 1; // gain is already corrected for STUB API! TODO:
				int offset = PANELSIZE * det_id;
				float *dest = img->buf + offset;

				for (int i = PANELSIZE; i > 0; i--) {
					*(dest++) = *(detData++) * gain;
				}
			}
		}

		cur_tags[det_id] += TAG_INCREMENT;
//		printf("cur_tags[%d] = %d\n", det_id, cur_tags[det_id]);
		nanosleep(&LONG_WAIT, NULL);
	}

	free(pDataStBuf);
	free(pWorkBuf);
	
	printf("Child%d: Finished\n", det_id);
}

int main(int argc, const char * argv[])
{	
   	printf("Cheetah for SACLA Online API version 221011\n\n");
	printf("Takanori Nakane, 2014-2022\n");
	printf("based on the work by Anton Barty, 21 January 2014\n");
	
	// Input data file and Cheetah configuration file
	char	filename[1024];
	char	cheetahini[1024];

	if (argc != 3 && argc != 4) {
		printf("Usage: %s input.h5 setting.ini [tag_increment; 1 or 2]\n", argv[0]);
		return -1;
	}

	// Take configuration from command line arguments
	strcpy(filename,argv[1]);
	strcpy(cheetahini,argv[2]);

	if (argc == 4) {
		TAG_INCREMENT = atoi(argv[3]);
		if (TAG_INCREMENT != 2 && TAG_INCREMENT != 1) {
			printf("tag_increment must be 1 or 2.\n");
			return -1;
		}
	}
    
	printf("Program name: %s\n",argv[0]);
	printf("Input data file: %s\n", filename);
	printf("Cheetah .ini file: %s\n", cheetahini);
	
	/*
	 *	Initialise Cheetah
	 */
	printf("Setting up Cheetah...\n");
	static uint32_t ntriggers = 0;
	static long frameNumber = 0;
	static cGlobal cheetahGlobal;
	static time_t startT = 0;
	time(&startT);
	strcpy(cheetahGlobal.configFile, cheetahini);
	strncpy(cheetahGlobal.cxiFilename, "output.h5", MAX_FILENAME_LENGTH);
	cheetahInit(&cheetahGlobal);
	cheetahGlobal.runNumber = -1;    

	/*
	 *	Initialize API
	 */
//	ol_initialize_dummy(filename);
    
	std::vector<std::string> detIDList, detIDListAll;

	struct da_string_array *det_ids;
	int n_det;
        da_alloc_string_array(&det_ids);
 	ol_read_detid_list(det_ids);
	da_getsize_string_array(&n_det, det_ids);
	for (int j = 0; j < n_det; j++) {
		char *val;
		da_getstring_string_array(&val, det_ids, j);
		detIDListAll.push_back(val);
		free(val);		
	}
	da_destroy_string_array(&det_ids);

	printf("API: ol_readDetIDList returned %d detectors.\n", (int)detIDListAll.size());
	for (int j = 0; j < detIDListAll.size(); j++) {
		printf("detID %d is %s\n", j, detIDListAll[j].c_str());
	}

	// Since this is online API, we assume there are no 'reconst' detectors.
	for (int j = 0; j < detIDListAll.size(); j++) {
		if (strncmp(detIDListAll[j].c_str(), "MPCCD-8", 7) == 0) {
			detIDList.push_back(detIDListAll[j]);
		}
	}
	   
	int ndet = detIDList.size();
	printf("Found %d MPCCD panels.\n", ndet);
	
	int sockIDs[NDET] = {};
	for (int detID = 0; detID < ndet; detID++) {
		ol_connect(&sockIDs[detID], detIDList[detID].c_str());
		printf("API: ol_connect for det %s (id %d) returned sockID %d\n", detIDList[detID].c_str(), detID, sockIDs[detID]);
	}

	int datasize, worksize;
	char *pDataStBufs[NDET] = {}, *pWorkBufs[NDET] = {};
	for (int detID = 0; SINGLE_THREAD && detID < ndet; detID++) {
		ol_get_data_size(&datasize, &worksize, sockIDs[detID]);
		printf("API: ol_getDataSize for %d returned dataSize = %d, workSize = %d\n", detID, datasize, worksize);

		pDataStBufs[detID] = (char*)malloc(datasize);
		pWorkBufs[detID] = (char*)malloc(worksize);
	}

	/*
	 *	Initialize threads
	 */
	
	pthread_t threads[NDET];
	pthread_mutex_t mutex_map;
	std::map <int, Image_Data*> images;
	int cur_tags[NDET];
	pthread_mutex_init(&mutex_map, NULL);

	for (int i = 0; i < ndet; i++) {
		cur_tags[i] = -1;
	}

	if (!SINGLE_THREAD) {
		for (int i = 0; i < ndet; i++) {
			Thread_Data *td = (Thread_Data*)malloc(sizeof(Thread_Data));
			td->detID = i;
			td->socketID = sockIDs[i];
			td->mutex_map = &mutex_map;
			td->images = &images;
			td->cur_tags = cur_tags;
			td->photon_energy = cheetahGlobal.defaultPhotonEnergyeV;
			pthread_create(threads + i, NULL, thread, td);
		}
	}

    /*
     * Create a buffer for holding the detector image data from all 8 panels
     */
    long    fs_one = XSIZE;
    long    ss_one = YSIZE;
    long    fs = fs_one;
    long    ss = 8 * ss_one;
    long    nn = fs*ss;
    float  *buffer = (float*)malloc(nn * sizeof(float));
    hsize_t dims[2];
    dims[0] = ss;
    dims[1] = fs;
    
	bool failed = false;
	int tagID = -1;
    int runNumber;

	while (tagID < 0 && !SINGLE_THREAD) {
		tagID = cur_tags[PRIMARY_DET];// TODO? - TAG_INCREMENT; 
		printf("MainThread: Waiting for the first image.\n");
		nanosleep(&SHORT_WAIT, NULL);
	}
	printf("MainThread: First image to process = %d\n", tagID);

	while (!failed) {
		bool image_ready = false;
		bool skip = false;
		Image_Data *img = NULL;
		std::map<int, Image_Data*>::const_iterator it;
		int wanted_tag = tagID + TAG_INCREMENT;

		if (!SINGLE_THREAD) {
//			 printf("MainThread: Waiting image %d\n", wanted_tag);
			
			// wait till an image become ready
			while (!image_ready) {
				image_ready = true;
				for (int i = 0; i < detIDList.size(); i++) {
					if (cur_tags[i] <= wanted_tag) {
						nanosleep(&SHORT_WAIT, NULL);
						image_ready = false;
						break;
					}
				}
			}
			printf("MainThread: Image %d is ready. Delay = %d frame(s).\n", wanted_tag, cur_tags[PRIMARY_DET] - wanted_tag);

			pthread_mutex_lock(&mutex_map);
			it = images.find(wanted_tag);
			if (it == images.end()) {
				printf("ERROR: Couldn't get image data for tag %d. This should not happen!\n", wanted_tag);
				skip = true;
			} else {
				img = it->second;

				if (img->error) {
					printf("ERROR: Image %d has its error flag set. skipped.\n", wanted_tag);
					skip = true;
				}
				images.erase(wanted_tag);
				pthread_mutex_unlock(&mutex_map);

				if (cur_tags[PRIMARY_DET] > wanted_tag + 3000) { // delay of 100 sec
					printf("ERROR: Image %d came too late. skipped.\n", wanted_tag);
					skip = true;
				}
			}

			if (!skip) {
				memcpy(buffer, img->buf, sizeof(float) * IMGSIZE);
				runNumber = img->run;
//				printf("MainThread: Removed image %d\n", wanted_tag);
			}

			if (img != NULL) {
				free(img);
			}
			tagID = wanted_tag;

		} else { // SINGLE THREAD
			for (int detID = 0; detID < ndet; detID++) {
				if (tagID == -1) wanted_tag = -1; // first image
				int err = ol_collect_det_data(pDataStBufs[detID], pWorkBufs[detID], &tagID, sockIDs[detID], wanted_tag, datasize, worksize);
				if (err < 0) {
					printf("ERROR: Couldn't get image data for tag %d.\n", wanted_tag);
					failed = true;
					break;
				}
				float *detData;
				ol_read_det_data(&detData, pDataStBufs[detID], 0);

				float gain;
				ol_read_abs_gain(&gain, pDataStBufs[detID], 0);
				gain *= 3.65 / 0.1 / cheetahGlobal.defaultPhotonEnergyeV;
		
				int offset = PANELSIZE * detID;
				for (int i = 0; i < PANELSIZE; i++) {
					buffer[offset + i] = detData[i] * gain;
				}
			}

			ol_read_run_num(&runNumber, pDataStBufs[0], 0);
		}
		
		if (skip) continue;
		if (failed) break;
		
		printf("Processing event: run = %d tag = %d energy = %f eV\n", runNumber, tagID, cheetahGlobal.defaultPhotonEnergyeV);
		frameNumber++;

		if (runNumber != cheetahGlobal.runNumber) {
			cheetahGlobal.runNumber = runNumber;
			cheetahNewRun(&cheetahGlobal);
			cheetahGlobal.writeInitialLog(true);
		}
		
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
            
		/*
		 *	Cheetah: Create a new eventData structure in which to place all information
		 */
		cEventData	*eventData;
		eventData = cheetahNewEvent(&cheetahGlobal);
		ntriggers++;
            
		/*
		 *  Cheetah: Populate event structure with meta-data
		 */
		eventData->frameNumber = tagID;
		eventData->runNumber = runNumber;
		eventData->nPeaks = 0;
		eventData->pumpLaserCode = 0;
		eventData->pumpLaserDelay = 0;
		eventData->photonEnergyeV = cheetahGlobal.defaultPhotonEnergyeV; // in eV
		eventData->wavelengthA = 12398 / eventData->photonEnergyeV; // 4.1357E-15 * 2.9979E8 * 1E10 / eV (A)
		eventData->pGlobal = &cheetahGlobal;
		eventData->fiducial = tagID; // must be unique
            
		/*  TODO: We can directly put into data_raw
		 *
		 *  Cheetah: Copy image data into
		 *  Raw data is currently hard-coded as UINT16_t, SACLA provides as float, so we have to loose precision :-(
		 */
		long    detID = 0;
		long    pix_nn = cheetahGlobal.detector[detID].pix_nn;
           
		for(long ii=0; ii<pix_nn; ii++) {
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
            
		/*
		 *	Cheetah: Process this event
		 */
		cheetahProcessEventMultithreaded(&cheetahGlobal, eventData);

		if (!SINGLE_THREAD) {
			// nanosleep(&SHORT_WAIT, NULL);
		}
	}
    	
	/*
	 *	Cheetah: Cleanly exit by closing all files, releasing memory, etc.
	 */
	cheetahExit(&cheetahGlobal);
 
    // actually, we should call join or use detached threads

	free(buffer);
	for (int detID = 0; SINGLE_THREAD && detID < ndet; detID++) {
		free(pDataStBufs[detID]);
		free(pWorkBufs[detID]);
	}
	
	time_t endT;
	time(&endT);
	double dif = difftime(endT,startT);
	std::cout << "time taken: " << dif << " seconds\n";
	
	// Clean exit
	std::cout << "Clean exit\n";
    return 0;
}

