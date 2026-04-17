#ifndef WRITER_THREAD_H
#define WRITER_THREAD_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "writer.h"
#include "options.h"
#include <atomic>
#include "hwy/contrib/thread_pool/futex.h"
#include <mutex>
#include <libdeflate.h>
#include "singleproducersingleconsumerlist.h"

using namespace std;

static constexpr int OFFSET_RING_SIZE = 512;

struct alignas(64) OffsetSlot {
    std::atomic<size_t> cumulative_offset{0};
    std::atomic<uint32_t> published_seq{UINT32_MAX};
};

class WriterThread{
public:
    WriterThread(Options* opt, string filename, bool isSTDOUT = false);
    ~WriterThread();

    void initWriter(string filename1, bool isSTDOUT = false);
    void initBufferLists();

    void cleanup();

    bool isCompleted();
    void output();
    void input(int tid, string* data);
    bool setInputCompleted();

    uint32_t bufferLength() {return mBufferLength;};
    void waitForBufferBelow(uint32_t limit) {
        for(;;) {
            uint32_t cur = mBufferLength.load(std::memory_order_acquire);
            if(cur <= limit) break;
            hwy::BlockUntilDifferent(cur, mBufferLength);
        }
    }
    string getFilename() {return mFilename;}
    bool isPwriteMode() {return mPwriteMode;}

private:
    void deleteWriter();
    void inputPwrite(int tid, string* data);
    void setInputCompletedPwrite();

private:
    Writer* mWriter1;
    Options* mOptions;
    string mFilename;

    bool mInputCompleted;
    std::atomic<uint32_t> mBufferLength;
    std::atomic<uint32_t> mWriterNotify;  // incremented to wake writer thread
    SingleProducerSingleConsumerList<string*>** mBufferLists;
    int mWorkingBufferList;

    // pwrite mode: parallel libdeflate gz compression + direct file write
    bool mPwriteMode;
    int mFd;
    OffsetSlot* mOffsetRing;
    std::atomic<size_t>* mNextSeq;
    libdeflate_compressor** mCompressors;
    char** mCompBufs;       // per-worker pre-allocated compress output buffers
    size_t* mCompBufSizes;  // per-worker buffer sizes
};

#endif
