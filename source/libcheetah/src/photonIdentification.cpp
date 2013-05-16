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


void calculateSignificanceMap(cEventData *eventData, cGlobal *global) {

  for (int detID = 0; detID < global->nDetectors; detID++) {
    uint16_t  *mask = eventData->detector[detID].pixelmask;
    float     *data = eventData->detector[detID].corrected_data;
    double    *sigma = global->detector[detID].darkSigmaMap;
    float *significanceMap = eventData->detector[detID].significanceMap;
    float *photonMap = eventData->detector[detID].photonMap;
    //float *significanceMap = global->detector[detID].significanceMap;
    //float *photonMap = global->detector[detID].photonMap;
    long *cumPhotonMap = global->detector[detID].cumPhotonMap;
    long      pix_nn = global->detector[detID].pix_nn;
    // Combine pixel options for pixels to be ignored
    uint16_t  pixel_options = PIXEL_IS_IN_PEAKMASK | PIXEL_IS_OUT_OF_RESOLUTION_LIMITS | PIXEL_IS_HOT | PIXEL_IS_BAD | PIXEL_IS_SATURATED | PIXEL_IS_MISSING;
    
    for (int i = 0; i < pix_nn; i++) {
      significanceMap[i] = fabs(data[i]) / fabs(sigma[i]);
      //printf("sig[%d] = %g = %g / %g \n", i, significanceMap[i], data[i], sigma[i]);
    }

    for (int i = 0; i < pix_nn; i++) {
      if (significanceMap[i] > global->sigPhotonThreshold) {
	// 30 ADU per photon
	photonMap[i] = data[i] / 30.;

#ifdef __GNUC__
	__sync_fetch_and_add(&(cumPhotonMap[i]),(long)photonMap[i]);
#else
	int lockThreads = global->detector[detID].cumPhoton_mutex;
	if (lockThreads) {pthread_mutex_lock(&global->gbbuffer_mutex);}
	global->detector[detID].cumPhoton_mutex[i] += (long)photonMap[i];
	if (lockThreads) {pthread_mutex_unlock(&global->cumPhoton_mutex);}
#endif
	//cumPhotonMap[i] += photonMap[i];
      } else {
	photonMap[i] = 0.;
      }
    }
  }
}


void calculatePhotonMap(cEventData *eventData, cGlobal *global) {

  for (int detID = 0; detID < global->nDetectors; detID++) {
    uint16_t  *mask = eventData->detector[detID].pixelmask;
    float     *data = eventData->detector[detID].corrected_data;
    float    *thresholdMap = global->detector[detID].thresholdMap;
    float *photonMap = eventData->detector[detID].photonMap;
    //float *significanceMap = global->detector[detID].significanceMap;
    //float *photonMap = global->detector[detID].photonMap;
    long *cumPhotonMap = global->detector[detID].cumPhotonMap;
    long      pix_nn = global->detector[detID].pix_nn;
    // Combine pixel options for pixels to be ignored
    uint16_t  pixel_options = PIXEL_IS_IN_PEAKMASK | PIXEL_IS_OUT_OF_RESOLUTION_LIMITS | PIXEL_IS_HOT | PIXEL_IS_BAD | PIXEL_IS_SATURATED | PIXEL_IS_MISSING;
    
    for (int i = 0; i < pix_nn; i++) {
      if (data[i] > thresholdMap[i]) {
	// 30 ADU per photon
	photonMap[i] = data[i] / 30.;

#ifdef __GNUC__
	__sync_fetch_and_add(&(cumPhotonMap[i]),(long)photonMap[i]);
#else
	int lockThreads = global->detector[detID].cumPhoton_mutex;
	if (lockThreads) {pthread_mutex_lock(&global->gbbuffer_mutex);}
	global->detector[detID].cumPhoton_mutex[i] += (long)photonMap[i];
	if (lockThreads) {pthread_mutex_unlock(&global->cumPhoton_mutex);}
#endif
	//cumPhotonMap[i] += photonMap[i];
      } else {
	photonMap[i] = 0.;
      }
    }
  }
}
