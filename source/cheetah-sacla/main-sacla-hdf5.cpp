//
//  main.cpp
//  cheetah-sacla-hdf5
//
//  Originally created by Anton Barty on 20/1/14.
//  Copyright (c) 2014 Anton Barty. All rights reserved.
//
//  Modified by Takanori Nakane to correct relative gains
//  and work in 0.1 photon units
//
//  gaincal file was not used for the sake of CPU cache efficiency.

#include <iostream>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "cheetah.h"
#include "sacla-hdf5-reader.h"

int main(int argc, const char * argv[])
{	
   	printf("Simple SACLA HDF5 file parser\n");
	printf("Anton Barty, 21 January 2014\n");
	printf("Takanori Nakane, 2014-2015\n");
	
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
	 *	Open SACLA HDF5 file
	 *	Read file header and information
	 */
    SACLA_h5_info_t SACLA_header = {};
    SACLA_HDF5_ReadHeader(filename, &SACLA_header);
    
    /*
     * Create a buffer for holding the detector image data from all 8 panels
     */
    long    fs_one = 512;
    long    ss_one = 1024;
    long    nn_one = fs_one*ss_one;
    long    fs = fs_one;
    long    ss = 8*ss_one;
    long    nn = fs*ss;
    float   *buffer = (float*) calloc(nn, sizeof(float));
    hsize_t dims[2];
    dims[0] = ss;
    dims[1] = fs;
    
    // Loop through runs in this HDF5 file
    for(long runID=0; runID<SACLA_header.nruns; runID++) {
        printf("Processing run: %s\n", SACLA_header.run_string[runID]);
        runNumber = SACLA_header.run_number[runID];
        
        // Gather detector fields and event tags for this run
        SACLA_HDF5_Read2dDetectorFields(&SACLA_header, runID);
        SACLA_HDF5_ReadRunInfo(&SACLA_header, runID);
        
        // Loop through all events found in this run
        for(long eventID=200; eventID<SACLA_header.nevents; eventID++) { // DEBUG!
            int tagID = atoi(SACLA_header.event_name[eventID] + 4); // "tag_######"
            printf("Processing event: tag = %d energy = %f eV\n", tagID, SACLA_header.actual_photon_energy_in_eV[eventID]);
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
			 *	SACLA: Read next image
			 */
            if (SACLA_HDF5_ReadImageRaw(&SACLA_header, runID, eventID, buffer, nn_one) < 0) {
				continue;
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
            eventData->photonEnergyeV = SACLA_header.actual_photon_energy_in_eV[eventID];        // in eV
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
            
			
            /*
             *  SACLA: Dump buffer to a trivial HDF5 file (to see whether the file reader is working at all)
             */
            /*
			  char    outfile[1024];
			  hid_t   outfile_id;
			  sprintf(outfile,"/data/scratch/sacla/%s.h5", SACLA_header.event_name[eventID]);
			  printf("Writing to temporary file: %s\n",outfile);
			  outfile_id = H5Fcreate(outfile,  H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
			  H5LTmake_dataset_float(outfile_id, "data", 2, dims, buffer );
			  H5Fclose(outfile_id);
			*/
            
        }
    }
    
	
	// Clean up stale IDs and exit
    SACLA_HDF5_cleanup(&SACLA_header);
    
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

