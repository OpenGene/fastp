// in-memory pooled Reads to reduce new/delete operations
// for each thread, one SISC list is used

#ifndef READ_POOL_H
#define READ_POOL_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include "options.h"
#include "singleproducersingleconsumerlist.h"

using namespace std;

class ReadPool{
public:
    ReadPool(Options* opt);
    ~ReadPool();
    bool input(int tid, Read* r);
    Read* getOne();
    void initBufferLists();
    void cleanup();
    size_t size();

private:
	void updateFullStatus();

private:
    Options* mOptions;
    SingleProducerSingleConsumerList<Read*>** mBufferLists;
    size_t mProduced;
    size_t mConsumed;
    unsigned long mLimit;
    bool mIsFull;
};

#endif