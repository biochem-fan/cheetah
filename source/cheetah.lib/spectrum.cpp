//
//  spectrum.cpp
//  cheetah
//
//  Created by Richard Bean on 2/27/13.
//
//


#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <hdf5.h>
#include <stdlib.h>

#include "detectorObject.h"
#include "cheetahGlobal.h"
#include "cheetahEvent.h"
#include "cheetahmodules.h"
#include "data2d.h"


void integrateSpectrum(cEventData *eventData, cGlobal *global) {
	// proceed if event is a 'hit', spectrum data exists & spectrum required
	int hit = eventData->hit;
	int opalfail = eventData->specFail;
	int specWidth = eventData->specWidth;
	int specHeight = eventData->specHeight;
	
	int spectra = global->espectrum1D;
	
	if(global->generateDarkcal && !opalfail && spectra){
		eventData->energySpectrumExist = 1;
		genSpectrumBackground(eventData,global,specWidth,specHeight);
	}
	if(hit && !opalfail && spectra){
		eventData->energySpectrumExist = 1;
		integrateSpectrum(eventData,global,specWidth,specHeight);
	}
	return;
}


void integrateSpectrum(cEventData *eventData, cGlobal *global, int specWidth,int specHeight) {
	// integrate spectrum into single line and output to event data
	float PIE = 3.141;
	float ttilt = tanf(global->espectiltang*PIE/180);
	int opalindex;
	int newind;

	for (long i=0; i<specHeight; i++) {
		for (long j=0; j<specWidth; j++) {
			newind = i + (int) ceilf(j*ttilt);        // index of the integrated array, must be integer,!
			if (newind >= 0 && newind < specHeight) {
				opalindex = i*specWidth + j;   // index of the 2D camera array
				eventData->energySpectrum1D[newind]+=eventData->specImage[opalindex];
				if (global->espectrumDarkSubtract) {
					eventData->energySpectrum1D[newind]-=global->espectrumDarkcal[opalindex];
				}
			}
		}
	}
	return;
}


void integrateRunSpectrum(cEventData *eventData, cGlobal *global) {
	
	// Update integrated run spectrum
	if(eventData->hit && eventData->energySpectrumExist) {
		pthread_mutex_lock(&global->espectrumRun_mutex);
		for (long i=0; i<global->espectrumLength; i++) {
			global->espectrumRun[i] += eventData->energySpectrum1D[i];
		}
		pthread_mutex_unlock(&global->espectrumRun_mutex);
	}

	// Update spectrum hit counter
	if(eventData->energySpectrumExist && !global->generateDarkcal) {
		pthread_mutex_lock(&global->nespechits_mutex);
		global->nespechits++;
		pthread_mutex_unlock(&global->nespechits_mutex);
	}
	return;
}


void genSpectrumBackground(cEventData *eventData, cGlobal *global, int specWidth, int specHeight) {
	// Generate background for spectrum detector
	int spectrumpix = specWidth*specHeight;

	pthread_mutex_lock(&global->espectrumBuffer_mutex);
	for (int i=0; i<spectrumpix; i++) {
		global->espectrumBuffer[i]+=eventData->specImage[i];
	}
	pthread_mutex_unlock(&global->espectrumBuffer_mutex);
	pthread_mutex_lock(&global->nespechits_mutex);
	global->nespechits++;
	pthread_mutex_unlock(&global->nespechits_mutex);
	return;
}


void saveIntegratedRunSpectrum(cGlobal *global) {

	int     spectrumpix = global->espectrumWidth*global->espectrumLength;
	double *espectrumDark = (double*) calloc(spectrumpix, sizeof(double));
	double *espectrumScale = (double*) calloc(global->espectrumLength, sizeof(double));
	char	filename[1024];
	int     maxindex = 0;
	int     evspread = global->espectrumSpreadeV;
	double  pixincrement = (double) evspread/global->espectrumLength;
	double  beamAveV = global->meanPhotonEnergyeV;
	double  eVoffset;

	// compute spectrum camera darkcal and save to HDF5
	if(global->generateDarkcal){
		pthread_mutex_lock(&global->espectrumRun_mutex);
		pthread_mutex_lock(&global->nespechits_mutex);
		for(int i=0; i<spectrumpix; i++) {
			espectrumDark[i] = global->espectrumBuffer[i]/global->nespechits;
		}

		sprintf(filename,"r%04u-energySpectrum-darkcal.h5", global->runNumber);
		printf("Saving energy spectrum darkcal to file: %s\n", filename);
        
		writeSimpleHDF5(filename, espectrumDark, global->espectrumWidth, global->espectrumLength, H5T_NATIVE_DOUBLE);

		free(espectrumDark);
		return;
	}

	// find maximum of run integrated spectum array and save both to HDF5
	pthread_mutex_lock(&global->espectrumRun_mutex);
	pthread_mutex_lock(&global->nespechits_mutex);

	for (int i=0; i<global->espectrumLength; i++) {
		if (global->espectrumRun[i] > global->espectrumRun[maxindex]) {
			maxindex = i;
		}
	}
	eVoffset = beamAveV - maxindex*pixincrement;
	for (int i=0; i<global->espectrumLength; i++) {
		espectrumScale[i]=i*pixincrement + eVoffset;
	}

	sprintf(filename,"r%04u-integratedEnergySpectrum.h5", global->runNumber);
	printf("Saving run-integrated energy spectrum: %s\n", filename);

	writeSpectrumInfoHDF5(filename, espectrumScale, global->espectrumRun, global->espectrumLength, H5T_NATIVE_DOUBLE, &maxindex, 1, H5T_NATIVE_INT);

	pthread_mutex_unlock(&global->espectrumRun_mutex);
	pthread_mutex_unlock(&global->nespechits_mutex);
	return;
}

void readSpectrumDarkcal(cGlobal *global, char *filename) {

	int spectrumpix = global->espectrumLength*global->espectrumWidth;

	// Do we need a darkcal file?
	if (global->espectrumDarkSubtract == 0){
		return;
	}
	
	// Check if a darkcal file has been specified
	if ( strcmp(filename,"") == 0 ){
		printf("spectrum camera Darkcal file path was not specified.\n");
		exit(1);
	}
	
	// Check whether file exists!
	FILE* fp = fopen(filename, "r");
	if (fp) 	// file exists
		fclose(fp);
	else {		// file doesn't exist
		printf("\tspecified energy spectrum Darkcal file does not exist: %s\n",filename);
		printf("\tAborting...\n");
		exit(1);
	}

	printf("Reading energy spectrum Darkcal file:\n");
	printf("\t%s\n",filename);
	
	// Read darkcal data from file
	cData2d		temp2d;
	temp2d.readHDF5(filename);
	// Copy into darkcal array
	for(long i=0; i<spectrumpix; i++) {
		global->espectrumDarkcal[i] =  temp2d.data[i];
	}
	return;
}

void readSpectrumEnergyScale(cGlobal *global, char *filename) {
	
	double*     energyscale = (double *) calloc(global->espectrumLength, sizeof(double));
	char        groupname[1024];
	char        fieldname[1024];
	hid_t       file_id;
	hid_t       datagroup_id;
	hid_t       dataset_id;
	hid_t       dataspace_id;
	hid_t       datatype_id;
	H5T_class_t dataclass;
	size_t      size;
	
	int ndims;
	
	sprintf(groupname, "energySpectrum");
	sprintf(fieldname, "runIntegratedEnergyScale");
	// Check if an energy scale calibration file has been specified
	if ( strcmp(filename,"") == 0 ){
		printf("spectrum energy scale calibration file path was not specified\n");
		printf("spectra will be output with default (0) energy scale\n");
		return;
	}
	
	// Check whether file exists!
	FILE* fp = fopen(filename, "r");
	if (fp) 	// file exists
		fclose(fp);
	else {		// file doesn't exist
		printf("specified energy scale calibration file does not exist: %s\n",filename);
		printf("spectra will be output with default (0) energy scale\n");
		return;
	}
	
	printf("Reading energy spectrum scale calibration file:\n");
	printf("\t%s\n",filename);
	
	// Open the file
	file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
	if(file_id < 0){
		printf("ERROR: Could not open file %s\n",filename);
		printf("spectra will be output with default (0) energy scale\n");
		return;
	}
	
	// Open the dataset
	datagroup_id = H5Gopen1(file_id, groupname);
	dataset_id = H5Dopen1(datagroup_id, fieldname);
	dataspace_id = H5Dget_space(dataset_id);
	
	// Test if correct dimensions / size
	ndims = H5Sget_simple_extent_ndims(dataspace_id);
	if(ndims != 1) {
		printf("the specified file does not have the correct dimensions for energy scale calibration, ndims=%i\n",ndims);
		printf("spectra will be output with default (0) energy scale\n");
		return;
	}
	hsize_t dims[ndims];
	H5Sget_simple_extent_dims(dataspace_id,dims,NULL);
	if (!dims[0]==1 || !dims[1]==global->espectrumLength) {
		printf("the specified file does not have the correct dimensions for energy scale calibration\n");
		printf("spectra will be output with default (0) energy scale\n");
		return;
	}
	
	datatype_id =  H5Dget_type(dataset_id);
	dataclass = H5Tget_class(datatype_id);
	size = H5Tget_size(datatype_id);
		
	H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, energyscale);
	for(int i=0; i<global->espectrumLength; i++) {
		global->espectrumScale[i] = energyscale[i];
	}
	free(energyscale);
	
	// Close and cleanup
	H5Dclose(dataset_id);
	H5Gclose(datagroup_id);
	
	// Cleanup stale IDs
	hid_t ids[256];
	int n_ids = H5Fget_obj_ids(file_id, H5F_OBJ_ALL, 256, ids);
	for (long i=0; i<n_ids; i++ ) {
		
		hid_t id;
		H5I_type_t type;
		id = ids[i];
		type = H5Iget_type(id);
		if ( type == H5I_GROUP )
			H5Gclose(id);
		if ( type == H5I_DATASET )
			H5Dclose(id);
		if ( type == H5I_DATASPACE )
			H5Sclose(id);
		//if ( type == H5I_DATATYPE )
		//	H5Dclose(id);
	}
	
	H5Fclose(file_id);
	printf("energy spectrum scale calibration file read successfull:\n");
	return;
}