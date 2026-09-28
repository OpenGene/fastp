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
#include <condition_variable>
#include <chrono>
#include <libdeflate.h>
#include "singleproducersingleconsumerlist.h"

using namespace std;

static constexpr int OFFSET_RING_SIZE = 512;

struct alignas(64) OffsetSlot {
    std::atomic<size_t> cumulative_offset{0};
    std::atomic<size_t> published_seq{SIZE_MAX};
};

// Holds one pending output chunk, keyed by its pack's global sequence number
// rather than by which worker produced it -- lets packs be distributed to
// workers out of strict round-robin order (see assignQueueForRound) while
// output is still reassembled in original read order.
struct alignas(64) OutputSlot {
    std::atomic<string*> data{nullptr};
    std::atomic<size_t> seq{SIZE_MAX};
};

class WriterThread{
public:
    WriterThread(Options* opt, string filename, bool isSTDOUT = false);
    ~WriterThread();

    void initWriter(string filename1, bool isSTDOUT = false);
    void initOutputRing();

    void cleanup();

    bool isCompleted();
    void output();
    // tid identifies the calling worker (stable for that worker's lifetime,
    // used only to index per-worker compression resources in pwrite mode);
    // seq is the pack's global sequence number and determines output order.
    void input(int tid, size_t seq, string* data);
    bool setInputCompleted();

    long bufferLength() {return mBufferLength;};
    string getFilename() {return mFilename;}
    bool isPwriteMode() {return mPwriteMode;}

private:
    void deleteWriter();
    void inputPwrite(int tid, size_t seq, string* data);
    void setInputCompletedPwrite();

private:
    Writer* mWriter1;
    Options* mOptions;
    string mFilename;

    bool mInputCompleted;
    atomic_long mBufferLength;
    OutputSlot* mOutputRing;
    size_t mOutputRingSize;
    std::atomic<size_t> mNextExpectedSeq;
    std::mutex mOutputMtx;
    std::condition_variable mOutputCV;

    // pwrite mode: parallel libdeflate gz compression + direct file write
    bool mPwriteMode;
    int mFd;
    OffsetSlot* mOffsetRing;
    std::atomic<size_t> mMaxPublishedSeq;
    libdeflate_compressor** mCompressors;
    char** mCompBufs;       // per-worker pre-allocated compress output buffers
    size_t* mCompBufSizes;  // per-worker buffer sizes
};

#endif
