#ifndef DUPLICATE_H
#define DUPLICATE_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include "options.h"
#include "common.h"
#include <atomic>

using namespace std;

class Duplicate{
public:
    Duplicate(Options* opt);
    ~Duplicate();

    bool checkRead(Read* r1);
    bool checkPair(Read* r1, Read* r2);
    void seq2intvector(const char* data, int len, uint64* output, int posOffset = 0);

    double getDupRate();

private:
    void initPrimeArrays();
    bool applyBloomFilter(uint64* positions);

private:
    Options* mOptions;
    uint64 mBufLenInBits;
    uint64 mBufLenInBytes;
    uint32 mBufNum;
    atomic_uchar* mDupBuf;
    uint64* mPrimeArrays;
    atomic_ulong mTotalReads;
    atomic_ulong mDupReads;
    uint64 mOffsetMask;
    
};

#endif