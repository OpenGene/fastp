#include "readpool.h"
#include "util.h"
#include <memory.h>
#include <unistd.h>
#include "common.h"

ReadPool::ReadPool(Options* opt){
    mOptions = opt;
    initBufferLists();
    mLimit = PACK_SIZE * PACK_IN_MEM_LIMIT;
    if(mOptions->interleavedInput)
        mLimit *= 2;
    mIsFull = false;
    mProduced = 0;
    mConsumed = 0;
}

ReadPool::~ReadPool() {
    cleanup();
}

bool ReadPool::input(int tid, Read* data) {
    if(mIsFull)
        return false;

    mBufferLists[tid]->produce(data);
    mProduced++;
    if((mProduced & 0xFF) == 0)
        updateFullStatus();
    return true;
}

void ReadPool::cleanup() {
    for(int t=0; t<mOptions->thread; t++) {
        while(mBufferLists[t]->canBeConsumed()) {
            Read* r = mBufferLists[t]->consume();
            mConsumed++;
            delete r;
        }
        delete mBufferLists[t];
    }
    delete[] mBufferLists;
}

void ReadPool::initBufferLists() {
    mBufferLists = new SingleProducerSingleConsumerList<Read*>*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++) {
        mBufferLists[t] = new SingleProducerSingleConsumerList<Read*>();
    }
}

void ReadPool::updateFullStatus() {
    mIsFull = size() > mLimit;
}

size_t ReadPool::size() {
    size_t total = 0;
    for(int t=0; t<mOptions->thread; t++) {
        total += mBufferLists[t]->size();
    }
    return total;
}

Read* ReadPool::getOne() {
    for(int t=0; t<mOptions->thread; t++) {
        if(mBufferLists[t]->canBeConsumed()) {
            Read* r = mBufferLists[t]->consume();
            mConsumed++;
            if((mConsumed & 0xFF) == 0)
                updateFullStatus();
            return r;
        }
    }
    return NULL;
}