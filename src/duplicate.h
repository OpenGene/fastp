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
    void addRecord(uint32 key, uint64 kmer32);

    // make histogram and get duplication rate
    static double statAll(vector<Duplicate*>& list, int* hist, int histSize);

private:
    Options* mOptions;
    int mKeyLenInBase;
    int mKeyLenInBit;
    uint64* mDups;
    uint16* mCounts;
    
};

#endif