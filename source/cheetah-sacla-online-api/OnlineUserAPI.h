/**
* @file OnlineUserAPI.h
* @par Copyright RIKEN
*/

#ifndef ONLINEUSERAPI_H_INCLUDE
#define ONLINEUSERAPI_H_INCLUDE
#include <string>
#include <vector>

// newest tag
#define OL_NEWESTTAGDATA  -1

/**
* error code (returned by function)
*/
enum OnlineAPI_ErrorCode{
	OL_ERR_TAGDATAGONE			= -10000,
	OL_ERR_WRONGDATASIZE		= -10010,
	OL_ERR_IMCOMPATIBLEDATASIZE	= -10011,
	OL_ERR_SMALLREADDATASIZE	= -10012,
	OL_ERR_SMALLDATABUFSIZE		= -10013,
	OL_ERR_SOCKET				= -10100,
	OL_ERR_SOCKOPT				= -10101,
	OL_ERR_GETHOSTBYNAME		= -10102,
	OL_ERR_CONNECT				= -10103,
	OL_ERR_NOSOCKET				= -10104,
	OL_ERR_NOTAGLO				= -10105,
	OL_ERR_TCPSEND				= -10106,
	OL_ERR_TCPRECV				= -10107,
	OL_ERR_SMALLWORKBUF			= -10108,
	OL_ERR_RECV					= -10109,
	OL_ERR_NULLPTR				= -10110,
	OL_ERR_WRONGDATAIDX			= -10111,
	OL_ERR_ARRYBUFALLOC			= -10200,
	OL_ERR_DATABUFALLOC			= -10201,
	OL_ERR_METABUFALLOC			= -10202,
	OL_ERR_WORKBUFALLOC			= -10203,
	OL_ERR_EVENTDATAALLOC		= -10204,
	OL_ERR_BUFALLOC				= -10205,
	OL_ERR_NOCALIBMATHOD		= -10300,
	OL_ERR_USERAPI_CHECK		= -500000
};

/**
* get version number of this API
*/
std::string ol_getAPIVersion();

/**
* read latest ID list of detector.
* @param[out] pDetIDList ID list of detector
*/
int ol_readDetIDList(std::vector<std::string> *pDetIDList);

/**
* connect to data handling server.
* @param[in] detID detector ID
* @param[out] pSockID socket ID to connect to data handling server
*/
int ol_connect(const char *detID, int *pSockID);

/**
* get size of data structure and working buffer.
* @param[in] sockID socket ID to connect to data handling server
* @param[out] pDataStSize buffer size for detector data structure
* @param[out] pWorkSize working buffer size
*/
int ol_getDataSize(int sockID, int *pDataStSize, int *pWorkSize);

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
int ol_collectDetData(int sockID, int tag, char *pDataStBuf, int dataStSize, char *pWorkBuf, int workSize, int *pTag);

/**
* read address of detector data.
* @param[in] pDataStBuf address of detector data structure
* @param[out] pData address of detector data
* @param[in] idx index of detector data (normally specify 0)
*/
int ol_readDetData(char *pDataStBuf, float **pData, int idx = 0);

/**
* read width and height of detector.
* @param[in] pDataStBuf address of detector data structure
* @param[out] pWidth width of detector data (pixel)
* @param[out] pHeight height of detector data (pixel)
* @param[in] idx index of detector data (normally specify 0)
*/
int ol_readDetSize(char *pDataStBuf, int *pWidth, int *pHeight, int idx = 0);

/**
* read run number.
* @param[in] pDataStBuf address of detector data structure
* @param[out] pNum run number
* @param[in] idx index of detector data (normally specify 0)
*/
int ol_readRunNum(char *pDataStBuf, int *pNum, int idx = 0);

/**
* read tag number.
* @param[in] pDataStBuf address of detector data structure
* @param[out] pNum tag number
* @param[in] idx index of detector data (normally specify 0)
*/
int ol_readTagNum(char *pDataStBuf, int *pNum, int idx = 0);

/**
* read absolute gain.
* @param[in] pDataStBuf address of detector data structure
* @param[out] pGain absolute gain
* @param[in] idx index of detector data (normally specify 0)
*/
int ol_readAbsGain(char *pDataStBuf, float *pGain, int idx = 0);


#endif // ONLINEUSERAPI_H_INCLUDE
