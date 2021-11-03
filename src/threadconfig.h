#ifndef THREAD_CONFIG_H
#define THREAD_CONFIG_H

#include "filterresult.h"
#include "options.h"
#include "singleproducersingleconsumerlist.h"
#include "stats.h"
#include "writer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

class ThreadConfig {
public:
  ThreadConfig(Options *opt, int threadId, bool paired = false);
  ~ThreadConfig();
  inline Stats *getPreStats1() { return mPreStats1; }
  inline Stats *getPostStats1() { return mPostStats1; }
  inline Stats *getPreStats2() { return mPreStats2; }
  inline Stats *getPostStats2() { return mPostStats2; }
  inline Writer *getWriter1() { return mWriter1; }
  inline Writer *getWriter2() { return mWriter2; }
  inline FilterResult *getFilterResult() { return mFilterResult; }

  void initWriter(string filename1);
  void initWriter(string filename1, string filename2);

  void addFilterResult(int result, int readNum);
  void addMergedPairs(int pairs);

  int getThreadId() { return mThreadId; }
  // for splitting output
  // increase mCurrentSplitReads by readNum, and check it with
  // options->split.size;
  void markProcessed(long readNum);
  void initWriterForSplit();
  bool canBeStopped();
  void cleanup();

  // input list
  void setInputList(SingleProducerSingleConsumerList<ReadPack *> *list);
  void setInputListPair(SingleProducerSingleConsumerList<ReadPack *> *left,
                        SingleProducerSingleConsumerList<ReadPack *> *right);
  SingleProducerSingleConsumerList<ReadPack *> *getLeftInput() {
    return mLeftInputList;
  }
  SingleProducerSingleConsumerList<ReadPack *> *getRightInput() {
    return mRightInputList;
  }

private:
  void deleteWriter();
  void writeEmptyFilesForSplitting();

private:
  Stats *mPreStats1;
  Stats *mPostStats1;
  Stats *mPreStats2;
  Stats *mPostStats2;
  Writer *mWriter1;
  Writer *mWriter2;
  Options *mOptions;
  FilterResult *mFilterResult;
  SingleProducerSingleConsumerList<ReadPack *> *mLeftInputList;
  SingleProducerSingleConsumerList<ReadPack *> *mRightInputList;

  // for spliting output
  int mThreadId;
  int mWorkingSplit;
  long mCurrentSplitReads;
  bool mCanBeStopped;
};

#endif