#ifndef SE_PROCESSOR_H
#define SE_PROCESSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include <cstdlib>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>
#include "options.h"
#include "threadconfig.h"
#include "filter.h"
#include "umiprocessor.h"
#include "writerthread.h"
#include "duplicate.h"
#include "singleproducersingleconsumerlist.h"
#include "readpool.h"

using namespace std;

typedef struct ReadRepository ReadRepository;

class SingleEndProcessor{
public:
    SingleEndProcessor(Options* opt);
    ~SingleEndProcessor();
    bool process();

private:
    bool processSingleEnd(ReadPack* pack, ThreadConfig* config);
    void readerTask();
    void processorTask(ThreadConfig* config);
    void initConfig(ThreadConfig* config);
    void initOutput();
    void closeOutput();
    void writerTask(WriterThread* config);
    void recycleToPool(int tid, Read* r);
    int pickLeastFullQueue();
    void releaseQueueSlot(int queueIndex);

private:
    Options* mOptions;
    atomic_bool mReaderFinished;
    alignas(128) atomic_int mFinishedThreads;
    Filter* mFilter;
    UmiProcessor* mUmiProcessor;
    WriterThread* mLeftWriter;
    WriterThread* mFailedWriter;
    Duplicate* mDuplicate;
    SingleProducerSingleConsumerList<ReadPack*>** mInputLists;
    size_t mPackReadCounter;
    alignas(128) atomic_long mPackProcessedCounter;
    long mPackInMemLimit;
    int mPackSize;
    ReadPool* mReadPool;
    std::mutex mBackpressureMtx;
    std::condition_variable mBackpressureCV;

    // Least-full-queue pack distribution (see #723 follow-up). SE has a
    // single reader thread producing to all worker queues, so unlike PE's
    // non-interleaved path there's no cross-thread coordination needed --
    // just pick the queue with the fewest in-flight packs each round.
    std::mutex mQueueAssignMtx;
    std::vector<int> mQueueDepth;
};


#endif