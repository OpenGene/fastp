#ifndef COMMON_H
#define COMMON_H

#define FASTP_VER "0.0.1"

typedef long int64;
typedef unsigned long uint64;

typedef int int32;
typedef unsigned int uint32;

typedef short int16;
typedef unsigned short uint16;

typedef char int8;
typedef unsigned char uint8;

#pragma pack(2) 


#pragma pack() 

// the limit of the queue to store the packs
// error may happen if it generates more packs than this number
static const int PACK_NUM_LIMIT  = 1000000;

// how many reads one pack has
static const int PACK_SIZE = 1000;

// if one pack is produced, but not consumed, it will be kept in the memory
// this number limit the number of in memory packs
// if the number of in memory packs is full, the producer thread should sleep
static const int PACK_IN_MEM_LIMIT = 100;

// if read number is more than this, warn it
static const int WARN_STANDALONE_READ_LIMIT = 10000;

// different filtering results
static const int FILTER_RESULT_TYPES = 6;
static const int PASS_FILTER = 0;
static const int FAIL_QUALITY = 1;
static const int FAIL_N_BASE = 2;
static const int FAIL_LENGTH = 3;
static const int FAIL_OVERLAP = 4;
static const int FAIL_POLY_X = 5;


#endif /* COMMON_H */
