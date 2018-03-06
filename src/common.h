#ifndef COMMON_H
#define COMMON_H

#define FASTP_VER "0.12.5"

#define _DEBUG false

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
static const int PACK_IN_MEM_LIMIT = 500;

// if read number is more than this, warn it
static const int WARN_STANDALONE_READ_LIMIT = 10000;

// different filtering results, bigger number means worse
// if r1 and r2 are both failed, then the bigger one of the two results will be recorded
// we reserve some gaps for future types to be added
static const int PASS_FILTER = 0;
static const int FAIL_POLY_X = 4;
static const int FAIL_OVERLAP = 8;
static const int FAIL_N_BASE = 12;
static const int FAIL_LENGTH = 16;
static const int FAIL_QUALITY = 20;

// how many types in total we support
static const int FILTER_RESULT_TYPES = 32;


#endif /* COMMON_H */
