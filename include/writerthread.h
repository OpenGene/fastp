#ifndef WRITER_THREAD_H
#define WRITER_THREAD_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "writer.h"
#include "options.h"
#include <atomic>
#include <mutex>

using namespace std;

class WriterThread{
public:
    WriterThread(Options* opt, string filename);
    ~WriterThread();

    void initWriter(string filename1);
    void initWriter(ofstream* stream);
    void initWriter(gzFile gzfile);

    void cleanup();

    bool isCompleted();
    void output();
    void input(char* data, size_t size);
    bool setInputCompleted();

    long bufferLength();
    string getFilename() {return mFilename;}

private:
    void deleteWriter();

private:
    Writer* mWriter1;
    Options* mOptions;
    string mFilename;

    // for spliting output
    bool mInputCompleted;
    atomic_long mInputCounter;
    atomic_long mOutputCounter;
    char** mRingBuffer;
    size_t* mRingBufferSizes;

    mutex mtx;

};

#endif