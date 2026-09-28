#ifndef PE_PROCESSOR_H
#define PE_PROCESSOR_H

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
#include "overlapanalysis.h"
#include "writerthread.h"
#include "duplicate.h"
#include "readpool.h"


using namespace std;

typedef struct ReadPairRepository ReadPairRepository;

class PairEndProcessor{
public:
    PairEndProcessor(Options* opt);
    ~PairEndProcessor();
    bool process();

private:
    bool processPairEnd(ReadPack* leftPack, ReadPack* rightPack, ThreadConfig* config);
    void readerTask(bool isLeft);
    void interleavedReaderTask();
    void processorTask(ThreadConfig* config);
    void initConfig(ThreadConfig* config);
    void initOutput();
    void closeOutput();
    void statInsertSize(Read* r1, Read* r2, OverlapResult& ov, int frontTrimmed1 = 0, int frontTrimmed2 = 0);
    int getPeakInsertSize();
    void writerTask(WriterThread* config);
    void recycleToPool1(int tid, Read* r);
    void recycleToPool2(int tid, Read* r);
    int pickLeastFullQueue();
    int assignQueueForRound(size_t round);
    void releaseQueueSlot(int queueIndex);

private:
    atomic_bool mLeftReaderFinished;
    atomic_bool mRightReaderFinished;
    alignas(128) atomic_int mFinishedThreads;
    Options* mOptions;
    Filter* mFilter;
    UmiProcessor* mUmiProcessor;
    atomic_long* mInsertSizeHist;
    WriterThread* mLeftWriter;
    WriterThread* mRightWriter;
    WriterThread* mUnpairedLeftWriter;
    WriterThread* mUnpairedRightWriter;
    WriterThread* mMergedWriter;
    WriterThread* mFailedWriter;
    WriterThread* mOverlappedWriter;
    Duplicate* mDuplicate;
    SingleProducerSingleConsumerList<ReadPack*>** mLeftInputLists;
    SingleProducerSingleConsumerList<ReadPack*>** mRightInputLists;
    size_t mLeftPackReadCounter;
    size_t mRightPackReadCounter;
    alignas(128) atomic_long mPackProcessedCounter;
    long mPackInMemLimit;
    int mPackSize;
    ReadPool* mLeftReadPool;
    ReadPool* mRightReadPool;
    atomic_bool shouldStopReading;
    std::mutex mBackpressureMtx;
    std::condition_variable mBackpressureCV;

    // Least-full-queue pack distribution: a worker's queue index is fixed at
    // construction, but which worker receives the *next* pack no longer has
    // to be strict round-robin. mQueueDepth tracks an approximate in-flight
    // pack count per worker, guarded by mQueueAssignMtx since it's
    // read-modify-written by up to two reader threads at once.
    std::mutex mQueueAssignMtx;
    std::vector<int> mQueueDepth;
    // The non-interleaved path reads left/right from two independent threads,
    // but a worker's mLeftInputLists[t]/mRightInputLists[t] pair must receive
    // matching pack #k on both sides. mQueueAssignments is a small ring cache
    // so whichever side reaches round k first picks the queue and the other
    // side reuses that exact choice, instead of each side picking separately
    // and risking a left/right mismatch. Sized comfortably above the maximum
    // possible left/right divergence, which mPackInMemLimit bounds.
    struct QueueAssignmentSlot { long round = -1; int queueIndex = 0; };
    std::vector<QueueAssignmentSlot> mQueueAssignments;
    size_t mQueueAssignRingSize;
};


#endif