#include "duplicate.h"
#include "overlapanalysis.h"
#include <memory.h>
#include <math.h>

const int PRIME_ARRAY_LEN = 300;

Duplicate::Duplicate(Options* opt) {
    mOptions = opt;

    // 1G mem required
    mBufLenInBytes = 1L <<29;
    mBufNum = 2;

    // if deduplication is enabled, we double the buffer size and buffer number to make more accurate deduplication
    // this will take 4x memory (4G)
    if(mOptions->duplicate.dedup) {
        mBufLenInBytes *= 2;
        mBufNum *= 2;
    }
    mBufLenInBits = mBufLenInBytes << 3;
    mDupBuf = new uint8[mBufLenInBytes * mBufNum];
    memset(mDupBuf, 0, sizeof(uint8) * mBufLenInBytes * mBufNum);

    mPrimeArrays = new uint64[mBufNum * PRIME_ARRAY_LEN];
    memset(mPrimeArrays, 0, sizeof(uint64) * mBufNum * PRIME_ARRAY_LEN);
    initPrimeArrays();

    mTotalReads = 0;
    mDupReads = 0;
}

void Duplicate::initPrimeArrays() {
    uint64 number = 1000;
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
            number += 100000;
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
                base = 11;
                break;
            case 'T':
                base = 31;
                break;
            case 'C':
                base = 71;
                break;
            case 'G':
                base = 101;
                break;
            default:
                base = 131;
        }
        for(int i=0; i<mBufNum; i++) {
            int offset = (p+posOffset)*mBufNum + i;
            offset %= (mBufNum*PRIME_ARRAY_LEN);
            output[i] += mPrimeArrays[offset] * (base + (p+posOffset));
        }
    }
}

bool Duplicate::checkRead(Read* r) {
    uint64* positions = new uint64[mBufNum];

    // init
    for(int i=0; i<mBufNum; i++)
        positions[i] = 0;
    int len = r->length();
    seq2intvector(r->mSeq.mStr.c_str(), len, positions);
    bool isDup = applyBloomFilter(positions);
    delete[] positions;

    mTotalReads++;
    if(isDup)
        mDupReads++;

    return isDup;
}

bool Duplicate::checkPair(Read* r1, Read* r2) {
    uint64* positions = new uint64[mBufNum];
    
    // init
    for(int i=0; i<mBufNum; i++)
        positions[i] = 0;
    seq2intvector(r1->mSeq.mStr.c_str(), r1->length(), positions);
    seq2intvector(r2->mSeq.mStr.c_str(), r2->length(), positions, r1->length());
    bool isDup = applyBloomFilter(positions);
    delete[] positions;

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

        isDup = isDup && (mDupBuf[i * mBufLenInBytes + bytePos] & byte);
        mDupBuf[i * mBufLenInBytes + bytePos] |= byte;
    }
    return isDup;
}

double Duplicate::getDupRate() {
    if(mTotalReads == 0)
        return 0.0;
    return (double)mDupReads/(double)mTotalReads;
}