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

#include "cheetah.h"
#include "OnlineUserAPI.h"

void ol_initialize_dummy(char* filename);

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
	
	int sockIDs[8] = {};
	for (int detID = 0; detID < 8; detID++) {
		ol_connect(detIDList[detID].c_str(), &sockIDs[detID]);
		printf("API: ol_connect for det %s returned sockID %d\n", detIDList[detID].c_str(), sockIDs[detID]);
	}

	int datasize, worksize;
	ol_getDataSize(sockIDs[0], &datasize, &worksize);
	printf("API: ol_getDataSize returned dataSize = %d, workSize = %d\n", datasize, worksize);

	char *pDataStBufs[8] = {}, *pWorkBufs[8] = {};
	for (int detID = 0; detID < 8; detID++) {
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
    float  *buffer = (float*) calloc(nn, sizeof(float));
    hsize_t dims[2];
    dims[0] = ss;
    dims[1] = fs;
    
	int runID = 0;

	for (int tag_loop = 0; tag_loop < 500; tag_loop++) {
		int tagID;

		for (int detID = 0; detID < 8; detID++) {
			ol_collectDetData(sockIDs[detID], -1, pDataStBufs[detID], datasize, pWorkBufs[detID], worksize, &tagID);
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
	
	time_t endT;
	time(&endT);
	double dif = difftime(endT,startT);
	std::cout << "time taken: " << dif << " seconds\n";
	
	// Clean exit
	std::cout << "Clean exit\n";
    return 0;
}

