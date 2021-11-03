#ifndef SE_PROCESSOR_H
#define SE_PROCESSOR_H

#include "duplicate.h"
#include "filter.h"
#include "options.h"
#include "read.h"
#include "readpool.h"
#include "singleproducersingleconsumerlist.h"
#include "threadconfig.h"
#include "umiprocessor.h"
#include "writerthread.h"
#include <condition_variable>
#include <cstdlib>
#include <mutex>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <thread>

using namespace std;

typedef struct ReadRepository ReadRepository;

class SingleEndProcessor {
public:
  SingleEndProcessor(Options *opt);
  ~SingleEndProcessor();
  bool process();

private:
  bool processSingleEnd(ReadPack *pack, ThreadConfig *config);
  void readerTask();
  void processorTask(ThreadConfig *config);
  void initConfig(ThreadConfig *config);
  void initOutput();
  void closeOutput();
  void writerTask(WriterThread *config);
  void recycleToPool(int tid, Read *r);

private:
  Options *mOptions;
  atomic_bool mReaderFinished;
  atomic_int mFinishedThreads;
  Filter *mFilter;
  UmiProcessor *mUmiProcessor;
  WriterThread *mLeftWriter;
  WriterThread *mFailedWriter;
  Duplicate *mDuplicate;
  SingleProducerSingleConsumerList<ReadPack *> **mInputLists;
  size_t mPackReadCounter;
  atomic_long mPackProcessedCounter;
  ReadPool *mReadPool;
};

#endif