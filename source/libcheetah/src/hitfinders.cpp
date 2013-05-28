#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <hdf5.h>
#include <stdlib.h>

#include "detectorObject.h"
#include "cheetahGlobal.h"
#include "cheetahEvent.h"
#include "median.h"
#include "hitfinders.h"
#include "peakfinders.h"


/*
 *	Various flavours of hitfinder
 *		0 - Everything is a hit
 *		1 - Number of pixels above ADC threshold
 *		2 - Total intensity above ADC threshold
 *		3 - Count Bragg peaks
 *		4 - Use TOF
 *		5 - Depreciated and no longer exists
 *		6 - Experimental - find peaks by SNR criteria
 *              7 - Laser on event code (usually EVR41)
 *              8 - Use TOF peak
 *              28  - Statistical hit finder
 *              29  - Reference based hitfinder
 *              30 - As 28 but with gmd
 *              31 - Threshold the photon count with a halo map
 */
int  hitfinder(cEventData *eventData, cGlobal *global){
	
	// Dereference stuff
	int	    detID = global->hitfinderDetector;
	int		hit=0;
	
	/*
	 *	Default values for some metrics
	 */
	eventData->peakNpix = 0;
	eventData->peakTotal = 0;
	eventData->peakResolution = 0;
	eventData->peakDensity = 0;
      
 	/*
	 *	Use one of various hitfinder algorithms
	 */
 	switch(global->hitfinderAlgorithm) {
		
	case 0 :	// Everything is a hit. Used for converting xtc to hdf
	  hit = 1;
	  break;	

	case 1 :	// Count the number of pixels above ADC threshold
	  hit = hitfinder1(global,eventData,detID);
	  break;	
	  
	case 2 :	// Integrated intensity above ADC threshold
	  hit = hitfinder2(global,eventData,detID);
	  break;
			
	case 3 : 	// Count number of Bragg peaks
	  hit = peakfinder3(global,eventData, detID);			
	  break;	

	case 4 :	// Use TOF signal to find hits
	  hit = hitfinder4(global,eventData,detID);
	  break;
						
	case 6 : 	// Count number of Bragg peaks
	  hit = peakfinder6(global,eventData, detID);
	  break;
            
	case 7 : 	// Return laser on event code
	  hit = eventData->laserEventCodeOn;
	  eventData->nPeaks = eventData->laserEventCodeOn;
	  break;

	case 8 :	// Use TOF signal, maximum peak, to find hits
	  hit = hitfinder8(global,eventData,detID);
	  break;

	case 9 :	// Use TOF signal, maximum peak, excluding classical htis
	  hit = hitfinder8(global,eventData,detID) && !(hitfinder1(global,eventData,detID));
	  break;

	case 28 :
	  hit = hitfinder28(global, eventData, detID);
	  break;

	case 29:
	  hit = hitfinder29(global, eventData, detID);
	  break;

	case 30:
	  hit = hitfinder30(global, eventData, detID);
	  break;

	case 31:
	  hit = hitfinder31(global, eventData, detID);
	  break;
			
	default :
	  printf("Unknown hit finding algorithm selected: %i\n", global->hitfinderAlgorithm);
	  printf("Stopping in confusion.\n");
	  exit(1);
	  break;
			
	}
	
	// Statistics on the peaks, for certain hitfinders
	if( eventData->nPeaks > 1 &&
	   ( global->hitfinderAlgorithm == 3 || global->hitfinderAlgorithm == 5 || global->hitfinderAlgorithm == 6 ) ) {
		   
		long	np;
		long  kk;
		float	resolution;
		float	resolutionA;	
		float	cutoff = 0.95;
		long		pix_nn = global->detector[detID].pix_nn;
		long		asic_nx = global->detector[detID].asic_nx;
		long		asic_ny = global->detector[detID].asic_ny;
		long	*inx = (long *) calloc(pix_nn, sizeof(long));
		long	*iny = (long *) calloc(pix_nn, sizeof(long));


		np = eventData->nPeaks;
		if(np >= global->hitfinderNpeaksMax) 
		   np = global->hitfinderNpeaksMax; 
		kk = (long) floor(cutoff*np);
	
	

		// Pixel radius resolution (bigger is better)
		float *buffer1 = (float*) calloc(global->hitfinderNpeaksMax, sizeof(float));
		for(long k=0; k<np; k++) 
			buffer1[k] = eventData->peak_com_r_assembled[k];
		resolution = kth_smallest(buffer1, np, kk);		   
		eventData->peakResolution = resolution;
		free(buffer1);
	
		// Resolution to real space (in Angstrom)
		// Crystallographic resolution d = lambda/sin(theta)
		float z = global->detector[0].detectorZ;
		float dx = global->detector[0].pixelSize;
		double r = sqrt(z*z+dx*dx*resolution*resolution);
		double sintheta = dx*resolution/r;
		resolutionA = eventData->wavelengthA/sintheta;
		eventData->peakResolutionA = resolutionA;

	
		if(resolution > 0) {
			float	area = (3.141*resolution*resolution)/(asic_ny*asic_nx);
			eventData->peakDensity = (cutoff*np)/area;
		}
		free(inx); 			
		free(iny);	

	   
	}
	
	// Update central hit counter
	if(hit) {
		pthread_mutex_lock(&global->nhits_mutex);
		global->nhits++;
		global->nrecenthits++;
		pthread_mutex_unlock(&global->nhits_mutex);
	}
	

	
	return(hit);
	
}

void integratePixAboveThreshold(float *data,uint16_t *mask,long pix_nn,float ADC_threshold,uint16_t pixel_options,long *nat,float *tat){

  *nat = 0;
  *tat = 0.0;

  for(long i=0;i<pix_nn;i++){
    if(isNoneOfBitOptionsSet(mask[i],pixel_options) && data[i] > ADC_threshold){
      *tat += data[i];
      *nat += 1;
    }
  }
}


/*
 *	Hitfinder #1
 *	Hit if number of pixels above ADC threshold exceeds MinPixCount
 */

int hitfinder1(cGlobal *global, cEventData *eventData, long detID){

  int       hit = 0;
  long      nat = 0;
  float     tat = 0.;
  uint16_t  *mask = eventData->detector[detID].pixelmask;
  float     *data = eventData->detector[detID].corrected_data;
  long	    pix_nn = global->detector[detID].pix_nn;  
  float     ADC_threshold = global->hitfinderADC;
  // Combine pixel options for pixels to be ignored
  uint16_t  pixel_options = PIXEL_IS_IN_PEAKMASK | PIXEL_IS_OUT_OF_RESOLUTION_LIMITS | PIXEL_IS_HOT | PIXEL_IS_BAD | PIXEL_IS_SATURATED | PIXEL_IS_MISSING | PIXEL_IS_IN_HALO;
  
  integratePixAboveThreshold(data,mask,pix_nn,ADC_threshold,pixel_options,&nat,&tat);
  eventData->peakTotal = tat;
  eventData->peakNpix = nat;
  eventData->nPeaks = nat;

  if(nat >= global->hitfinderMinPixCount){
    hit = 1;
  }
  return hit;
}


/*
 *	Hitfinder #2
 *	Hit if integrated value of pixels above ADC threshold exceeds TAT threshold
 */

int hitfinder2(cGlobal *global, cEventData *eventData, long detID){
  
  int       hit = 0;
  long      nat = 0;
  float     tat = 0.;
  uint16_t  *mask = eventData->detector[detID].pixelmask;
  float     *data = eventData->detector[detID].corrected_data;
  long	    pix_nn = global->detector[detID].pix_nn;  
  float     ADC_threshold = global->hitfinderADC;
  // Combine pixel options for pixels to be ignored
  uint16_t  pixel_options = PIXEL_IS_IN_PEAKMASK | PIXEL_IS_OUT_OF_RESOLUTION_LIMITS | PIXEL_IS_HOT | PIXEL_IS_BAD | PIXEL_IS_SATURATED | PIXEL_IS_MISSING;
  
  integratePixAboveThreshold(data,mask,pix_nn,ADC_threshold,pixel_options,&nat,&tat);
  eventData->peakTotal = tat;
  eventData->peakNpix = nat;
  eventData->nPeaks = nat;

  if(tat >= global->hitfinderTAT){
    hit = 1;
  }
  return hit;
}

int hitfinder4(cGlobal *global,cEventData *eventData,long detID){
  int hit = 0;
  long		pix_nn = global->detector[detID].pix_nn;
  uint16_t      *mask = eventData->detector[detID].pixelmask;
  long	nat = 0;
  long	counter;
  float	total;
  float	mingrad = global->hitfinderMinGradient*2;
  mingrad *= mingrad;
  
  nat = 0;
  counter = 0;
  total = 0.0;

  /*
   *	Create a buffer for image data so we don't nuke the main image by mistake 
   */
  float *temp = (float*) calloc(pix_nn, sizeof(float));
  memcpy(temp, eventData->detector[detID].corrected_data, pix_nn*sizeof(float));
	

  // combine pixelmask bits
  uint16_t combined_pixel_options = PIXEL_IS_IN_PEAKMASK | PIXEL_IS_OUT_OF_RESOLUTION_LIMITS | PIXEL_IS_HOT | PIXEL_IS_BAD | PIXEL_IS_SATURATED;
	
  /*
   *	Apply masks
   *	(multiply data by 0 to ignore regions)
   */
  for(long i=0;i<pix_nn;i++){
    temp[i] *= isNoneOfBitOptionsSet(mask[i], combined_pixel_options);
  }
		  
  if ((global->hitfinderUseTOF==1) && (eventData->TOFPresent==1)){
    double total_tof = 0.;
    for(int i=global->hitfinderTOFMinSample; i<global->hitfinderTOFMaxSample; i++){
      total_tof += eventData->TOFVoltage[i];
    }
    if (total_tof > global->hitfinderTOFThresh)
      hit = 1;
  }
  // Use cspad threshold if TOF is not present 
  else {
    for(long i=0;i<pix_nn;i++){
      if(temp[i] > global->hitfinderADC){
	nat++;
      }
    }
    if(nat >= global->hitfinderMinPixCount)
      hit = 1;
  }
  free(temp);
  return hit;
}



int hitfinder8(cGlobal *global,cEventData *eventData,long detID){
  int hit = 0;
  long		pix_nn = global->detector[detID].pix_nn;
  uint16_t      *mask = eventData->detector[detID].pixelmask;
  long	nat = 0;
  long	counter;
  float	total;
  float	mingrad = global->hitfinderMinGradient*2;
  mingrad *= mingrad;
  
  nat = 0;
  counter = 0;
  total = 0.0;

  /*
   *	Apply masks
   *	(multiply data by 0 to ignore regions)
   */

  if ((eventData->TOFPresent==1)){
    /*
    const int nback = 3;
    float olddata[nback];
    for (int k = 0; k < nback; k++)
      {
	olddata[k] = NAN;
      }
    for(int i=global->hitfinderTOFMinSample; i<global->hitfinderTOFMaxSample; i++){
      olddata[i % nback] = eventData->TOFVoltage[i];
      double sum = 0;
      for (int k = 0; k < nback; k++)
	{
	  sum += olddata[k];
	}
      if (sum < global->hitfinderTOFThresh * nback) hit = 1;
    }
    */
    float minSignal = 1.e-30;
    const int nback = 3;
    float olddata[nback];
    for (int k = 0; k < nback; k++)
      {
	olddata[k] = NAN;
      }
    for(int i=global->hitfinderTOFMinSample; i<global->hitfinderTOFMaxSample; i++){
      olddata[i % nback] = eventData->TOFVoltage[i];
      double sum = 0;
      for (int k = 0; k < nback; k++)
	{
	  sum += olddata[k];
	}
      if (sum < minSignal) minSignal=sum;
      //if (sum < global->hitfinderTOFThresh * nback) hit = 1;
    }
    if (minSignal < global->hitfinderTOFThresh * nback) hit = 1;
    eventData->tofIntegratedSignal = minSignal/(float)nback;
  } else {
    eventData->tofIntegratedSignal = 0.;
  }
  return hit;
}

int hitfinder28(cGlobal *global, cEventData *eventData, long detID) {
  
  int       hit = 0;
  uint16_t  *mask = eventData->detector[detID].pixelmask;
  float     *data = eventData->detector[detID].corrected_data;
  long	    pix_nn = global->detector[detID].pix_nn;  
  // Combine pixel options for pixels to be ignored
  uint16_t  pixel_options = PIXEL_IS_IN_PEAKMASK | PIXEL_IS_OUT_OF_RESOLUTION_LIMITS | PIXEL_IS_HOT | PIXEL_IS_BAD | PIXEL_IS_SATURATED | PIXEL_IS_MISSING;

  /*
  float photonSum = 0.;
  for (int i = 0; i<pix_nn; i++) {
    if (isNoneOfBitOptionsSet(mask[i],pixel_options) && (eventData->detector[detID].photonMap[i] > 0.)) {
      //photonSum += eventData->detector[detID].photonMap[i];
      photonSum += 1.;
    }
  }
  
  if (photonSum > global->totalPhotonsThreshold) {
    hit = 1;
  }
  */
  if (eventData->detector[detID].totalPhotons > global->totalPhotonsThreshold) {
    hit = 1;
  }
  return hit;
}


int hitfinder29(cGlobal *global, cEventData *eventData, long detID) {
  int hit = 0;
  uint16_t *mask = eventData->detector[detID].pixelmask;
  float *data = eventData->detector[detID].corrected_data;
  long pix_nn = global->detector[detID].pix_nn;
  uint16_t pixel_options = PIXEL_IS_IN_PEAKMASK | PIXEL_IS_OUT_OF_RESOLUTION_LIMITS | PIXEL_IS_HOT | PIXEL_IS_BAD | PIXEL_IS_SATURATED | PIXEL_IS_MISSING;

  float *photonMap = eventData->detector[detID].photonMap;
  float *sampleReference = global->detector[detID].sampleReference;
  float *backgroundReference = global->detector[detID].backgroundReference;

  float sampleProjection = 0.;
  float backgroundProjection = 0.;
  for (int i = 0; i < pix_nn; i++) {
    if (isNoneOfBitOptionsSet(mask[i],pixel_options) && (photonMap[i] > 0.)) {
      sampleProjection += sampleReference[i]*photonMap[i];
      backgroundProjection += backgroundReference[i]*photonMap[i];
    }
  }
  printf(">>> %g = %g / %g\n", sampleProjection / backgroundProjection, sampleProjection, backgroundProjection);
  if (sampleProjection / backgroundProjection > 1.) {
    hit = 1;
  }
  return hit;
}

int hitfinder30(cGlobal *global, cEventData *eventData, long detID) {
  int hit = 0;
  
  double gmd11 = eventData->gmd11;
  float totalPhotons = eventData->detector[detID].totalPhotons;
  
  // these numbers only make sense for original detector position
  float base = -8270.;
  float slope = 6459.;
  float sigma = 1110.;
  float threshold = 2.5;
  
  float expectedHaloSignal = base + slope*gmd11;
  if (gmd11 > 1. && (totalPhotons > expectedHaloSignal + sigma*threshold)) {
    hit = 1;
  }

  return hit;  
}

int hitfinder31(cGlobal *global, cEventData *eventData, long detID) {
  int hit = 0;
  int pix_nn = global->detector[detID].pix_nn;
  uint16_t *mask = eventData->detector[detID].pixelmask;
  uint16_t pixel_options = PIXEL_IS_IN_PEAKMASK | PIXEL_IS_OUT_OF_RESOLUTION_LIMITS | PIXEL_IS_HOT | PIXEL_IS_BAD | PIXEL_IS_SATURATED | PIXEL_IS_MISSING;
  float *photonMap = eventData->detector[detID].photonMap;
  float *haloMeanMap = global->detector[detID].haloMeanMap;
  float *haloSigmaMap = global->detector[detID].haloSigmaMap;
  float haloThreshold = global->haloThreshold;

  float sigma;
  float sigmaSum = 0.;
  for (int i = 0; i < pix_nn; i++) {
    if (isNoneOfBitOptionsSet(mask[i], pixel_options) && haloSigmaMap[i] > 0.){
      //printf("mean=%g, sigma=%g\n", haloMeanMap[i], haloSigmaMap[i]);
      sigma = (photonMap[i] - haloMeanMap[i]) / haloSigmaMap[i];
      /* do we add the sigmas or threshold on the sigmas? */
      /* try adding for now */
      /* try adding the squares to make it a log probability */
      sigmaSum+=pow(sigma,2);
    }
  }
  sigmaSum = sqrt(sigmaSum);
  eventData->detector[detID].haloSigma = sigmaSum;
  //printf("sigma sum = %g\n", sigmaSum);
  if (sigmaSum > haloThreshold) {
    hit = 1;
  }
  return hit;
}
