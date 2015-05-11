// Dummy implementation of SACLA online API for testing

#include "stdlib.h"

#include "OnlineUserAPI.h"
#include "../cheetah-sacla/sacla-hdf5-reader.h"

using namespace std;

typedef struct {
	int tag;
	int run;
	float gain;
	float buf[512 * 1024];
} DataStructure;

std::string ol_getAPIVersion() {
	string version = "OnlineAPI_DUMMY";
	return version;
}

SACLA_h5_info_t SACLA_header = {};
int SACLA_currentTag[8] = {};
float SACLA_buffer[8][512 * 8192] = {};

void ol_initialize_dummy(char* filename) {
    SACLA_HDF5_ReadHeader(filename, &SACLA_header);

    int runID = 0;
    int runNumber = SACLA_header.run_number[runID];
        
    // Gather detector fields and event tags for this run
    SACLA_HDF5_Read2dDetectorFields(&SACLA_header, runID);
    SACLA_HDF5_ReadRunInfo(&SACLA_header, runID);
}

/**
* read latest ID list of detector.
* @param[out] pDetIDList ID list of detector
*/
int ol_readDetIDList(std::vector<std::string> *pDetIDList) {
	char buf[10] = {};

	for (int i = 0; i < 8; i++) {
		buf[0] = '0' + i;
		pDetIDList->push_back(buf);
	}

	return 0;
}

/**
* connect to data handling server.
* @param[in] detID detector ID
* @param[out] pSockID socket ID to connect to data handling server
*/
int ol_connect(const char *detID, int *pSockID) {
	int tmp = detID[0] - '0';

	if (tmp < 0 || tmp >= 8) return -1;
	*pSockID = tmp;
	return 0;
}

/**
* get size of data structure and working buffer.
* @param[in] sockID socket ID to connect to data handling server
* @param[out] pDataStSize buffer size for detector data structure
* @param[out] pWorkSize working buffer size
*/
int ol_getDataSize(int sockID, int *pDataStSize, int *pWorkSize) {
	*pDataStSize = sizeof(DataStructure);
	*pWorkSize = sizeof(DataStructure);
	return 0;
}

/**
* collect detector data of specified tag.
* @param[in] sockID socket ID to connect to data handling server
* @param[in] tag tag number of detector data to acquire. latest tag : -1
* @param[out] pDataStBuf destination address of detector data structure
* @param[int] dataStSize buffer size for detector data structure
* @param[out] pWorkBuf destination address of working buffer
* @param[int] workSize working buffer size
* @param[out] pTag tag number acquired
*/
int ol_collectDetData(int sockID, int tag, char *pDataStBuf, int dataStSize, char *pWorkBuf, int workSize, int *pTag) {
    int runID = 0;
    
	if (sockID < 0 || sockID >= 8) return -1;
	if (dataStSize != sizeof(DataStructure)) return -1;

    if (tag == -1) {
		*pTag = SACLA_currentTag[sockID];
		SACLA_currentTag[sockID]++;
	} else {
		*pTag = tag;
	}
//	printf("API: ol_collectDetData sockID = %d, tag = %d, actual tag = %d\n", sockID, tag, *pTag);
	if (*pTag >= SACLA_header.nevents) return -1;

	if (*pTag >= 100 && false) { // DEBUG
		printf("API: on_collectDetData intentioanlly failing for DEBUG PURPOSE!\n");
		return -1;
	}

    SACLA_HDF5_ReadImageRaw(&SACLA_header, runID, *pTag, SACLA_buffer[sockID], 512 * 1024);
	int offset = 512 * 1024 * sockID;

	memcpy(((DataStructure*)pDataStBuf)->buf, SACLA_buffer[sockID] + offset, sizeof(float) * 512 * 1024);

	return 0;
}

/**
* read address of detector data.
* @param[in] pDataStBuf address of detector data structure
* @param[out] pData address of detector data
* @param[in] idx index of detector data (normally specify 0)
*/
int ol_readDetData(char *pDataStBuf, float **pData, int idx) {
	*pData = ((DataStructure*)pDataStBuf)->buf;

	return 0;
}

/**
* read width and height of detector.
* @param[in] pDataStBuf address of detector data structure
* @param[out] pWidth width of detector data (pixel)
* @param[out] pHeight height of detector data (pixel)
* @param[in] idx index of detector data (normally specify 0)
*/
int ol_readDetSize(char *pDataStBuf, int *pWidth, int *pHeight, int idx) {
	*pWidth = 512;
	*pHeight = 1024;

	return 0;
}

/**
* read run number.
* @param[in] pDataStBuf address of detector data structure
* @param[out] pNum run number
* @param[in] idx index of detector data (normally specify 0)
*/
int ol_readRunNum(char *pDataStBuf, int *pNum, int idx) {
	*pNum =  ((DataStructure*)pDataStBuf)->run;

	return 0;
}

/**
* read tag number.
* @param[in] pDataStBuf address of detector data structure
* @param[out] pNum tag number
* @param[in] idx index of detector data (normally specify 0)
*/
int ol_readTagNum(char *pDataStBuf, int *pNum, int idx) {
	*pNum = ((DataStructure*)pDataStBuf)->tag;

	return 0;
}

/**
* read absolute gain.
* @param[in] pDataStBuf address of detector data structure
* @param[out] pGain absolute gain
* @param[in] idx index of detector data (normally specify 0)
*/
int ol_readAbsGain(char *pDataStBuf, float *pGain, int idx) {
	*pGain = ((DataStructure*)pDataStBuf)->gain;

	return 0;
}
