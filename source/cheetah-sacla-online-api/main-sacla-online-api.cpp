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

#define NDET 8
#define IMGSIZE (512 * 1024 * 8)
#define SINGLE_THREAD false
#define TAG_INCREMENT 1 // 30 Hz

const struct timespec SHORT_WAIT = {0, 1E6}; // 1E6 = msec
const struct timespec LONG_WAIT = {0, 25E6};

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
	ol_getDataSize(socket_id, &datasize, &worksize);
	printf("API: ol_getDataSize returned dataSize = %d, workSize = %d\n", datasize, worksize);
	
	char *pDataStBuf = (char*)malloc(datasize);
	char *pWorkBuf = (char*)malloc(worksize);

	// Get the tag to start at
	
	if (det_id == 0) {
		int actual_tag = -1;
		int err = ol_collectDetData(socket_id, -1, pDataStBuf, datasize, pWorkBuf, worksize, &actual_tag);
		if (err < 0) {
			printf("ERROR: Failed to read the first image.\n");
			exit(-1);
		}

		while (actual_tag % TAG_INCREMENT != 0) {
			actual_tag++;
		}
		
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

		for (int ntry = 0; ntry < 30; ntry++) { // TODO: Check if this is enough
			pthread_mutex_lock(mutex_map);
			it = images->find(wanted_tag);
			if (it == images->end()) {
				if (det_id == 0) {
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
			int err = ol_collectDetData(socket_id, wanted_tag, pDataStBuf, datasize, pWorkBuf, worksize, &actual_tag);
			if (err < 0 || actual_tag != wanted_tag) {
				printf("Child%d: could not get tag %d. skipped. error code %d\n", det_id, wanted_tag, err);
				img->error = true;
				// This happens, for example, if we come too late (TAGDATAGONE -10000)
			} else {
				// printf("Child%d: Got image %d\n", det_id, tag);
				if (det_id == 0) {
					ol_readRunNum(pDataStBuf, &img->run);
				}

				float *detData;
				ol_readDetData(pDataStBuf, &detData);				
				int offset = 512 * 1024 * det_id;
				memcpy(img->buf + offset, detData, sizeof(float) * 512 * 1024);
			}
		}

		cur_tags[det_id] += TAG_INCREMENT;
		nanosleep(&LONG_WAIT, NULL);
	}

	free(pDataStBuf);
	free(pWorkBuf);
	
	printf("Child%d: Finished\n", det_id);
}

int main(int argc, const char * argv[])
{	
   	printf("Cheetah for SACLA Online API\n\n");
	printf("Takanori Nakane, 2014-2015\n");
	printf("based on the work by Anton Barty, 21 January 2014\n");
	
	// Input data file and Cheetah configuration file
	char	filename[1024];
	char	cheetahini[1024];

	if (argc != 3) {
		printf("Usage: %s input.h5 setting.ini\n", argv[0]);
		return -1;
	}

	// Take configuration from command line arguments
	strcpy(filename,argv[1]);
	strcpy(cheetahini,argv[2]);
    
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
    
    /*
	 *	Initialize API
	 */
	ol_initialize_dummy(filename);
    
    std::vector<std::string> detIDList;
	
	ol_readDetIDList(&detIDList);
	int ndet = detIDList.size();
	printf("API: ol_readDetIDList returned %d detectors.\n", ndet);
	
	int sockIDs[NDET] = {};
	for (int detID = 0; detID < ndet; detID++) {
		ol_connect(detIDList[detID].c_str(), &sockIDs[detID]);
		printf("API: ol_connect for det %s returned sockID %d\n", detIDList[detID].c_str(), sockIDs[detID]);
	}

	int datasize, worksize;
	ol_getDataSize(sockIDs[0], &datasize, &worksize);
	printf("API: ol_getDataSize returned dataSize = %d, workSize = %d\n", datasize, worksize);

	char *pDataStBufs[NDET] = {}, *pWorkBufs[NDET] = {};
	for (int detID = 0; SINGLE_THREAD && detID < ndet; detID++) {
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

	if (!SINGLE_THREAD) {
		for (int i = 0; i < ndet; i++) {
			cur_tags[i] = -1;

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
    long    fs_one = 512;
    long    ss_one = 1024;
    long    fs = fs_one;
    long    ss = 8*ss_one;
    long    nn = fs*ss;
    float  *buffer = (float*)malloc(nn * sizeof(float));
    hsize_t dims[2];
    dims[0] = ss;
    dims[1] = fs;
    
	bool failed = false;
	int tagID = -1; // TODO: how to initialize?
    int runNumber = 0;

	while (!failed) {
		bool image_ready = false;
		bool skip = false;
		Image_Data *img = NULL;
		std::map<int, Image_Data*>::const_iterator it;
		int wanted_tag = tagID + TAG_INCREMENT;

		if (!SINGLE_THREAD) {
			// printf("Main thread: Waiting image %d\n", wanted_tag);
			
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
			// printf("Main thread: Image %d is ready\n", wanted_tag);

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

				if (cur_tags[0] > wanted_tag + 600) { // delay of 20 sec
					printf("ERROR: Image %d came too late. skipped.\n", wanted_tag);
					skip = true;
				}
			}

			if (!skip) {
				memcpy(buffer, img->buf, sizeof(float) * 512 * 8192);
				runNumber = img->run;
				printf("Main thread: Removed image %d\n", wanted_tag);
			}

			if (img != NULL) {
				free(img);
			}

			// if (tagID > 1000) failed = true; // DEBUG
		} else { // SINGLE THREAD
			for (int detID = 0; detID < ndet; detID++) {
				int err = ol_collectDetData(sockIDs[detID], wanted_tag, pDataStBufs[detID], datasize, pWorkBufs[detID], worksize, &tagID);
				if (err < 0) {
					printf("ERROR: Couldn't get image data for tag %d.\n", wanted_tag);
					failed = true;
					break;
				}
				float *detData;
				ol_readDetData(pDataStBufs[detID], &detData);
//              ol_readAbsGain(pDataStBuf, &gain)

				int offset = 512 * 1024 * detID;				
				memcpy(buffer + offset, detData, sizeof(float) * 512 * 1024);				
			}

			ol_readRunNum(pDataStBufs[0], &runNumber);
		}
		tagID = wanted_tag;
		
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
			eventData->detector[detID].data_raw16[ii] = (uint16_t) lrint(buffer[ii]);
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

