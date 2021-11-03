#ifndef WRITER_THREAD_H
#define WRITER_THREAD_H

#include "options.h"
#include "singleproducersingleconsumerlist.h"
#include "writer.h"
#include <atomic>
#include <mutex>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

class WriterThread {
public:
  WriterThread(Options *opt, string filename);
  ~WriterThread();

  void initWriter(string filename1);
  void initBufferLists();

  void cleanup();

  bool isCompleted();
  void output();
  void input(int tid, string *data);
  bool setInputCompleted();

  long bufferLength() { return mBufferLength; };
  string getFilename() { return mFilename; }

private:
  void deleteWriter();

private:
  Writer *mWriter1;
  Options *mOptions;
  string mFilename;

  // for spliting output
  bool mInputCompleted;
  atomic_long mBufferLength;
  SingleProducerSingleConsumerList<string *> **mBufferLists;
  int mWorkingBufferList;
};

#endif