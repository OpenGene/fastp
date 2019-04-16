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
#include "options.h"
#include "threadconfig.h"
#include "filter.h"
#include "umiprocessor.h"
#include "writerthread.h"
#include "duplicate.h"

using namespace std;

struct ReadPack {
    Read** data;
    int count;
};

typedef struct ReadPack ReadPack;

struct ReadRepository {
    ReadPack** packBuffer;
    atomic_long readPos;
    atomic_long writePos;
    //std::mutex mtx;
    //std::mutex readCounterMtx;
    //std::condition_variable repoNotFull;
    //std::condition_variable repoNotEmpty;
};

typedef struct ReadRepository ReadRepository;

class SingleEndProcessor{
public:
    SingleEndProcessor(Options* opt);
    ~SingleEndProcessor();
    bool process();

private:
    bool processSingleEnd(ReadPack* pack, ThreadConfig* config);
    void initPackRepository();
    void destroyPackRepository();
    void producePack(ReadPack* pack);
    void consumePack(ThreadConfig* config);
    void producerTask();
    void consumerTask(ThreadConfig* config);
    void initConfig(ThreadConfig* config);
    void initOutput();
    void closeOutput();
    void writeTask(WriterThread* config);

private:
    Options* mOptions;
    ReadRepository mRepo;
    atomic_bool mProduceFinished;
    atomic_int mFinishedThreads;
    std::mutex mInputMtx;
    std::mutex mOutputMtx;
    Filter* mFilter;
    gzFile mZipFile;
    ofstream* mOutStream;
    UmiProcessor* mUmiProcessor;
    WriterThread* mLeftWriter;
    WriterThread* mFailedWriter;
    Duplicate* mDuplicate;
};


#endif