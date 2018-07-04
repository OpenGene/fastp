#include "duplicate.h"
#include "overlapanalysis.h"
#include <memory.h>
#include <math.h>

Duplicate::Duplicate(Options* opt) {
    mOptions = opt;
    mKeyLenInBase = mOptions->duplicate.keylen;
    mKeyLenInBit = 1<<(2*mKeyLenInBase);
    mDups = new uint64[mKeyLenInBit];
    memset(mDups, 0, sizeof(uint64)*mKeyLenInBit);
    mCounts = new uint16[mKeyLenInBit];
    memset(mCounts, 0, sizeof(uint16)*mKeyLenInBit);
    mGC = new uint8[mKeyLenInBit];
    memset(mGC, 0, sizeof(uint8)*mKeyLenInBit);
}

Duplicate::~Duplicate(){
    delete[] mDups;
    delete[] mCounts;
}

uint64 Duplicate::seq2int(const char* data, int start, int keylen, bool& valid) {
    uint64 ret = 0;
    for(int i=0; i<keylen; i++) {
        switch(data[start + i]) {
            case 'A':
                ret += 0;
                break;
            case 'T':
                ret += 1;
                break;
            case 'C':
                ret += 2;
                break;
            case 'G':
                ret += 3;
                break;
            default:
                valid = false;
                return 0;
        }
        // if it's not the last one, shift it by 2 bits
        if(i != keylen-1)
            ret <<= 2;
    }
    return ret;
}

void Duplicate::addRecord(uint32 key, uint64 kmer32, uint8 gc) {
    if(mCounts[key] == 0) {
        mCounts[key] = 1;
        mDups[key] = kmer32;
        mGC[key] = gc;
    } else {
        if(mDups[key] == kmer32)
            mCounts[key]++;
        else if(mDups[key] > kmer32) {
            mDups[key] = kmer32;
            mCounts[key] = 1;
            mGC[key] = gc;
        }
    }
}

void Duplicate::statRead(Read* r) {
    if(r->length() < 32)
        return;

    int start1 = 0;
    int start2 = max(0, r->length() - 32 - 5);

    const char* data = r->mSeq.mStr.c_str();
    bool valid = true;

    uint64 ret = seq2int(data, start1, mKeyLenInBase, valid);
    uint32 key = (uint32)ret;
    if(!valid)
        return;

    uint64 kmer32 = seq2int(data, start2, 32, valid);
    if(!valid)
        return;

    int gc = 0;

    // not calculated
    if(mCounts[key] == 0) {
        for(int i=0; i<r->length(); i++) {
            if(data[i] == 'C' || data[i] == 'T')
                gc++;
        }
    }

    gc = round(255.0 * (double) gc / (double) r->length());

    addRecord(key, kmer32, (uint8)gc);
}

void Duplicate::statPair(Read* r1, Read* r2) {
    if(r1->length() < 32 || r2->length() < 32)
        return;

    const char* data1 = r1->mSeq.mStr.c_str();
    const char* data2 = r2->mSeq.mStr.c_str();
    bool valid = true;

    uint64 ret = seq2int(data1, 0, mKeyLenInBase, valid);
    uint32 key = (uint32)ret;
    if(!valid)
        return;

    uint64 kmer32 = seq2int(data2, 0, 32, valid);
    if(!valid)
        return;

    int gc = 0;

    // not calculated
    if(mCounts[key] == 0) {
        for(int i=0; i<r1->length(); i++) {
            if(data1[i] == 'G' || data1[i] == 'C')
                gc++;
        }
        for(int i=0; i<r2->length(); i++) {
            if(data2[i] == 'G' || data2[i] == 'C')
                gc++;
        }
    }

    gc = round(255.0 * (double) gc / (double)( r1->length() + r2->length()));

    addRecord(key, kmer32, gc);
}

double Duplicate::statAll(int* hist, double* meanGC, int histSize) {
    long totalNum = 0;
    long dupNum = 0;
    int* gcStatNum = new int[histSize];
    memset(gcStatNum, 0, sizeof(int)*histSize);
    for(int key=0; key<mKeyLenInBit; key++) {
        int count = mCounts[key];
        double gc = mGC[key];

        if(count > 0) {
            totalNum += count;
            dupNum += count - 1;

            if(count >= histSize){
                hist[histSize-1]++;
                meanGC[histSize-1] += gc;
                gcStatNum[histSize-1]++;
            }
            else{
                hist[count]++;
                meanGC[count] += gc;
                gcStatNum[count]++;
            }
        }
    }

    for(int i=0; i<histSize; i++) {
        if(gcStatNum[i] > 0) {
            meanGC[i] = meanGC[i] / 255.0 / gcStatNum[i];
        }
    }

    delete[] gcStatNum;

    if(totalNum == 0)
        return 0.0;
    else
        return (double)dupNum / (double)totalNum;
}