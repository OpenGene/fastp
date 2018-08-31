#ifndef DUPLICATE_H
#define DUPLICATE_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include "options.h"
#include "common.h"

using namespace std;

class Duplicate{
public:
    Duplicate(Options* opt);
    ~Duplicate();

    void statRead(Read* r1);
    void statPair(Read* r1, Read* r2);
    uint64 seq2int(const char* data, int start, int keylen, bool& valid);
    void addRecord(uint32 key, uint64 kmer32, uint8 gc);

    // make histogram and get duplication rate
    double statAll(int* hist, double* meanGC, int histSize);

private:
    Options* mOptions;
    int mKeyLenInBase;
    int mKeyLenInBit;
    uint64* mDups;
    uint16* mCounts;
    uint8* mGC;
    
};

#endif