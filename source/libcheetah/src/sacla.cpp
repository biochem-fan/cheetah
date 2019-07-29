/*
 *  sacla.cpp
 *  SACLA specific additions to Cheetah
 *
 *  Created by Takanori Nakane, 2015-2019
 *  license: same as cheetah itself
 *
 */

/*
 *	Write out processed data to SACLA multi-event HDF5 format
 */

#include <stdio.h>
#include <math.h>
#include <hdf5.h>
#include <stdlib.h>

#include "detectorObject.h"
#include "cheetahGlobal.h"
#include "cheetahEvent.h"

void writeSACLA(cEventData *eventData, cGlobal *global) {
	// Update cleaned.txt (cf. saveFrame.cpp)
	pthread_mutex_lock(&global->framefp_mutex);
	fprintf(global->cleanedfp, "r%04u/%s/%s, %li, %i, %g, %g, %g, %g, %g\n",global->runNumber, eventData->eventSubdir, eventData->eventname, eventData->frameNumber, eventData->nPeaks, eventData->peakNpix, eventData->peakTotal, eventData->peakResolution, eventData->peakResolutionA, eventData->peakDensity);
	pthread_mutex_unlock(&global->framefp_mutex);

	const int detIndex = 0;
	// from saveFrame.cpp
	int16_t* corrected_data_int16 = (int16_t*)calloc(global->detector[detIndex].pix_nn, sizeof(int16_t));
			
	for (long i=0;i<global->detector[0].pix_nn;i++) {
		long tmp = lrint(eventData->detector[detIndex].data_detPhotCorr[i]);
		if (tmp < 0) {
			tmp = 0; 
		} else if (tmp > SHRT_MAX) {
			tmp = SHRT_MAX; 
		}
		corrected_data_int16[i] = (int16_t)tmp;
	}

	// This is a multi-event file, so a mutex is necessary
	pthread_mutex_lock(&global->saveCXI_mutex);
	hid_t file_id, dataset_id, dataspace_id, group_id;
	hsize_t dims[] = {global->detector[detIndex].pix_ny, global->detector[detIndex].pix_nx};

	file_id = H5Fopen(global->cxiFilename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0) { // create file for the first time. TODO: remove HDFlib warning
		file_id = H5Fcreate(global->cxiFilename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	}			

	char group_name[256];
	snprintf(group_name, 256, "/tag-%d", eventData->fiducial);
	group_id = H5Gcreate2(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	char dataset_name[256];
	dataspace_id = H5Screate(H5S_SCALAR);
	snprintf(dataset_name, 256, "%s/photon_energy_ev", group_name);
	dataset_id = H5Dcreate1(file_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &eventData->photonEnergyeV);
	H5Dclose(dataset_id);
	
	snprintf(dataset_name, 256, "%s/photon_wavelength_A", group_name);
	dataset_id = H5Dcreate1(file_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &eventData->wavelengthA);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);

	dataspace_id = H5Screate_simple(2, dims, NULL);
	snprintf(dataset_name, 256, "%s/data", group_name);
        
        hid_t dcpl = H5P_DEFAULT;
        if (global->h5compress != 0) {
                dcpl = H5Pcreate(H5P_DATASET_CREATE);
                //printf("global->h5compress = %d\n", global->h5compress);
                H5Pset_shuffle(dcpl);
                H5Pset_deflate(dcpl, global->h5compress);
                H5Pset_chunk(dcpl, 2, dims);
        }
        // was H5T_NATIVE_USHORT
	dataset_id = H5Dcreate2(file_id, dataset_name, H5T_STD_I16LE, dataspace_id, H5P_DEFAULT, dcpl, H5P_DEFAULT);

	H5Dwrite(dataset_id, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, corrected_data_int16);
	free(corrected_data_int16);

        H5Pclose(dcpl);
	H5Sclose(dataspace_id);
	H5Gclose(group_id);
	H5Dclose(dataset_id);
	H5Fclose(file_id);
	pthread_mutex_unlock(&global->saveCXI_mutex);
}
