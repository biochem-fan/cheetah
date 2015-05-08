//
//  cheetah-sacla-online-api
//
//  Originally created by Anton Barty on 20/1/14.
//  Copyright (c) 2014 Anton Barty. All rights reserved.
//
//  Modified by Takanori Nakane to read data from Online API
//
//  gaincal file was not used for the sake of CPU cache efficiency.

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

typedef struct {
  unsigned short buf[IMGSIZE];
  int error;
} Image_Data;

typedef struct {
  int detid;
} Thread_Data;

int acc = 0;
pthread_mutex_t mutex_map;
std::map <int, Image_Data*> images;
int cur_tags[NDET] = {}; // tag to look at NEXT
struct timespec small_wait = {0, 5E8}; // 1E6
struct timespec large_wait = {1, 0}; // 1E6

void* thread(void *thread_data) {
  int det_id = ((Thread_Data*)thread_data)->detid;
  free(thread_data);

  printf("Child thread %d: Created\n", det_id);

  while (1) {
    std::map<int, Image_Data*>::const_iterator it;
    Image_Data *img = NULL;
    int tag = cur_tags[det_id];

    while (1) {
      pthread_mutex_lock(&mutex_map);
      it = images.find(tag);
      if (it == images.end()) {
	if (det_id == 0) {
	  // Only thread 0 can allocate new buffer
	  img = (Image_Data*)calloc(1, sizeof(Image_Data));
	  images[tag] = img;
	  printf("Child thread %d: Allocated image %d\n", det_id, tag);
	}
      } else {
	img = it->second;
      }
      pthread_mutex_unlock(&mutex_map);

      if (img == NULL) {
	// wait till buffer is allocated
	nanosleep(&small_wait, NULL);
	printf("Child thread %d: Waiting\n", det_id);
      } else {
	break;
      }
    }
   
    // TODO: retrieve an image

    printf("Child thread %d: Got image %d\n", det_id, tag);
    cur_tags[det_id]++;
    nanosleep(&small_wait, NULL);
  }

  printf("Child thread %d: Finished\n", det_id);
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
    
    // Also for testing
    printf("Program name: %s\n",argv[0]);
    printf("Input data file: %s\n", filename);
    printf("Cheetah .ini file: %s\n", cheetahini);
	
	/*
	 *	Initialise Cheetah
	 */
	printf("Setting up Cheetah...\n");
	static uint32_t ntriggers = 0;
	static long frameNumber = 0;
    long runNumber = 0;
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
	printf("API: ol_readDetIDList returned %d detectors.\n", detIDList.size());
	
	int sockIDs[NDET] = {};
	for (int detID = 0; detID < detIDList.size(); detID++) {
		ol_connect(detIDList[detID].c_str(), &sockIDs[detID]);
		printf("API: ol_connect for det %s returned sockID %d\n", detIDList[detID].c_str(), sockIDs[detID]);
	}

	int datasize, worksize;
	ol_getDataSize(sockIDs[0], &datasize, &worksize);
	printf("API: ol_getDataSize returned dataSize = %d, workSize = %d\n", datasize, worksize);

	char *pDataStBufs[NDET] = {}, *pWorkBufs[NDET] = {};
	for (int detID = 0; detID < detIDList.size(); detID++) {
		pDataStBufs[detID] = (char*)malloc(datasize);
		pWorkBufs[detID] = (char*)malloc(worksize);
	}

    /*
     * Create a buffer for holding the detector image data from all 8 panels
     */
    long    fs_one = 512;
    long    ss_one = 1024;
    long    nn_one = fs_one*ss_one;
    long    fs = fs_one;
    long    ss = 8*ss_one;
    long    nn = fs*ss;
    float  *buffer = (float*)malloc(nn * sizeof(float));
    hsize_t dims[2];
    dims[0] = ss;
    dims[1] = fs;
    
	int runID = 0;
	bool failed = false;

	while (!failed) {
		int tagID;

		for (int detID = 0; detID < detIDList.size(); detID++) {
			int err = ol_collectDetData(sockIDs[detID], -1, pDataStBufs[detID], datasize, pWorkBufs[detID], worksize, &tagID);
			if (err < 0) {
				failed = true;
				break;
			}
			float *detData;
			ol_readDetData(pDataStBufs[detID], &detData);
//			ol_readRunNum(pDataStBuf, &run);
//          ol_readTagNum(pDataStBuf, &tag);
//          ol_readAbsGain(pDataStBuf, &gain)
			int offset = 512 * 1024 * detID;
			
			for (int i = 0; i < 512 * 1024; i++) {
				buffer[offset + i] = detData[i];
			}
		}

		if (failed) break;
		
		printf("Processing event: tag = %d energy = %f eV\n", tagID, 7000.0);
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
		eventData->photonEnergyeV = 7000;        // in eV
		eventData->wavelengthA = 12398 / eventData->photonEnergyeV; // 4.1357E-15 * 2.9979E8 * 1E10 / eV (A)
		eventData->pGlobal = &cheetahGlobal;
		eventData->fiducial = tagID; // must be unique
            
		/*
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
	}
    	
	/*
	 *	Cheetah: Cleanly exit by closing all files, releasing memory, etc.
	 */
	cheetahExit(&cheetahGlobal);

	free(buffer);
	for (int detID = 0; detID < detIDList.size(); detID++) {
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

