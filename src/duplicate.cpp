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

// Branchless base-to-hash-value lookup (A=7, T=222, C=74, G=31, else=13)
static const uint64 SEQ_HASH_VAL[256] = {
    13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, // 0x00-0x0F
    13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, // 0x10-0x1F
    13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, // 0x20-0x2F
    13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, // 0x30-0x3F
    13, 7,13,74,13,13,13,31,13,13,13,13,13,13,13,13, // 0x40-0x4F: A=7, C=74, G=31
    13,13,13,13,222,13,13,13,13,13,13,13,13,13,13,13, // 0x50-0x5F: T=222
    13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, // 0x60-0x6F
    13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, // 0x70-0x7F
    13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, // 0x80-0x8F
    13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, // 0x90-0x9F
    13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, // 0xA0-0xAF
    13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, // 0xB0-0xBF
    13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, // 0xC0-0xCF
    13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, // 0xD0-0xDF
    13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, // 0xE0-0xEF
    13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, // 0xF0-0xFF
};

void Duplicate::seq2intvector(const char* data, int len, uint64* output, int posOffset) {
    for(int p=0; p<len; p++) {
        uint64 base = SEQ_HASH_VAL[static_cast<unsigned char>(data[p])];
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