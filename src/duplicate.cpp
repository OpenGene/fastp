#include "duplicate.h"
#include "overlapanalysis.h"
#include <memory.h>
#include <math.h>
#include "util.h"

const int PRIME_ARRAY_LEN = 1<<9;

Duplicate::Duplicate(Options* opt) {
    mOptions = opt;

    // 1G mem required
    mBufLenInBytes = 1L <<29;
    mBufNum = 2;

    // the memory usage increases as the accuracy level increases
    // level 1: 1G
    // level 2: 2G
    // level 3: 4G
    // level 4: 8G
    // level 5: 16G
    // level 6: 32G
    switch(mOptions->duplicate.accuracyLevel) {
        case 1:
            break;
        case 2:
            mBufLenInBytes *= 2;
            break;
        case 3:
            mBufLenInBytes *= 2;
            mBufNum *= 2;
            break;
        case 4:
            mBufLenInBytes *= 4;
            mBufNum *= 2;
            break;
        case 5:
            mBufLenInBytes *= 8;
            mBufNum *= 2;
            break;
        case 6:
            mBufLenInBytes *= 8;
            mBufNum *= 4;
            break;
        default:
            break;
    }

    mOffsetMask = PRIME_ARRAY_LEN * mBufNum - 1;

    mBufLenInBits = mBufLenInBytes << 3;
    mDupBuf = new atomic_uchar[mBufLenInBytes * mBufNum];
    if(!mDupBuf) {
        error_exit("Out of memory, failed to allocate " + to_string(mBufLenInBytes * mBufNum) + " bytes buffer for duplication analysis, please reduce the dup_accuracy_level and try again.");
    }
    memset(mDupBuf, 0, sizeof(atomic_uchar) * mBufLenInBytes * mBufNum);

    mPrimeArrays = new uint64[mBufNum * PRIME_ARRAY_LEN];
    memset(mPrimeArrays, 0, sizeof(uint64) * mBufNum * PRIME_ARRAY_LEN);
    initPrimeArrays();

    mTotalReads = 0;
    mDupReads = 0;
}

void Duplicate::initPrimeArrays() {
    uint64 number = 10000;
    uint64 count = 0;
    while(count < mBufNum * PRIME_ARRAY_LEN) {
        number++;
        bool isPrime = true;
        for(uint64 i=2; i<=sqrt(number); i++) {
            if(number%i == 0) {
                isPrime = false;
                break;
            }
        }
        if(isPrime) {
            mPrimeArrays[count] = number;
            count++;
            number += 10000;
        }
    }
}

Duplicate::~Duplicate(){
    delete[] mDupBuf;
    delete[] mPrimeArrays;
}

void Duplicate::seq2intvector(const char* data, int len, uint64* output, int posOffset) {
    for(int p=0; p<len; p++) {
        uint64 base = 0;
        switch(data[p]) {
            case 'A':
                base = 7;
                break;
            case 'T':
                base = 222;
                break;
            case 'C':
                base = 74;
                break;
            case 'G':
                base = 31;
                break;
            default:
                base = 13;
        }
        for(int i=0; i<mBufNum; i++) {
            int offset = (p+posOffset)*mBufNum + i;
            offset &= mOffsetMask;
            output[i] += mPrimeArrays[offset] * (base + (p+posOffset));
        }
    }
}

bool Duplicate::checkRead(Read* r) {
    // max mBufNum is 8 (accuracy level 6), stack-allocate to avoid heap alloc per read
    uint64 positions[8] = {0};
    int len = r->length();
    seq2intvector(r->mSeq->c_str(), len, positions);
    bool isDup = applyBloomFilter(positions);

    mTotalReads++;
    if(isDup)
        mDupReads++;

    return isDup;
}

bool Duplicate::checkPair(Read* r1, Read* r2) {
    // max mBufNum is 8 (accuracy level 6), stack-allocate to avoid heap alloc per read
    uint64 positions[8] = {0};
    seq2intvector(r1->mSeq->c_str(), r1->length(), positions);
    seq2intvector(r2->mSeq->c_str(), r2->length(), positions, r1->length());
    bool isDup = applyBloomFilter(positions);

    mTotalReads++;
    if(isDup)
        mDupReads++;

    return isDup;
}

bool Duplicate::applyBloomFilter(uint64* positions) {
    bool isDup = true;
    for(int i=0; i<mBufNum; i++) {
        uint64 pos = positions[i] % mBufLenInBits;
        uint64 bytePos = pos >> 3;
        uint32 bitOffset = pos & 0x07;
        uint8 byte = (0x01) << bitOffset;

        //isDup = isDup && (mDupBuf[i * mBufLenInBytes + bytePos] & byte);
        uint8 ret = atomic_fetch_or(mDupBuf + i * mBufLenInBytes + bytePos, byte);
        isDup &= (ret & byte) != 0;
    }
    return isDup;
}

double Duplicate::getDupRate() {
    if(mTotalReads == 0)
        return 0.0;
    return (double)mDupReads/(double)mTotalReads;
}