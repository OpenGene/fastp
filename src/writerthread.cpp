#include "writerthread.h"
#include "util.h"
#include <memory.h>
#include <unistd.h>
#include <fcntl.h>
#include <cerrno>
#include <cstring>
#include <thread>
#include <chrono>

WriterThread::WriterThread(Options* opt, string filename, bool isSTDOUT){
    mOptions = opt;
    mWriter1 = NULL;
    mInputCompleted = false;
    mFilename = filename;

    mPwriteMode = !isSTDOUT && ends_with(filename, ".gz") && mOptions->thread > 1;
    mFd = -1;
    mOffsetRing = NULL;
    mNextSeq = NULL;
    mCompressors = NULL;
    mCompBufs = NULL;
    mCompBufSizes = NULL;
    mBufferLists = NULL;

    if (mPwriteMode) {
        mFd = open(mFilename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (mFd < 0)
            error_exit("Failed to open for pwrite: " + mFilename);
        mOffsetRing = new OffsetSlot[OFFSET_RING_SIZE];
        mNextSeq = new std::atomic<size_t>[mOptions->thread];
        for (int t = 0; t < mOptions->thread; t++)
            mNextSeq[t].store(t, std::memory_order_relaxed);
        mCompressors = new libdeflate_compressor*[mOptions->thread];
        for (int t = 0; t < mOptions->thread; t++)
            mCompressors[t] = libdeflate_alloc_compressor(mOptions->compression);
        size_t initBufSize = PACK_SIZE * 500;
        mCompBufs = new char*[mOptions->thread];
        mCompBufSizes = new size_t[mOptions->thread];
        for (int t = 0; t < mOptions->thread; t++) {
            mCompBufs[t] = new char[initBufSize];
            mCompBufSizes[t] = initBufSize;
        }
        mWorkingBufferList = 0;
        mBufferLength = 0;
        mWriterNotify = 0;
    } else {
        initWriter(filename, isSTDOUT);
        initBufferLists();
        mWorkingBufferList = 0;
        mBufferLength = 0;
        mWriterNotify = 0;
    }
}

WriterThread::~WriterThread() {
    cleanup();
}

bool WriterThread::isCompleted()
{
    if (mPwriteMode) return true;  // no writer thread needed
    return mInputCompleted && (mBufferLength==0);
}

bool WriterThread::setInputCompleted() {
    if (mPwriteMode) {
        setInputCompletedPwrite();
        mInputCompleted = true;
        return true;
    }
    mInputCompleted = true;
    for(int t=0; t<mOptions->thread; t++) {
        mBufferLists[t]->setProducerFinished();
    }
    // Wake the writer thread blocked in output() so it re-checks isCompleted()
    mWriterNotify.fetch_add(1, std::memory_order_release);
    hwy::WakeAll(mWriterNotify);
    return true;
}

void WriterThread::setInputCompletedPwrite() {
    // Acquire fence: synchronize with the release stores in inputPwrite()
    // so that all mNextSeq[t] writes from worker threads are visible here.
    std::atomic_thread_fence(std::memory_order_acquire);
    int W = mOptions->thread;
    size_t lastSeq = 0;
    bool anyProcessed = false;
    for (int t = 0; t < W; t++) {
        if (mNextSeq[t].load(std::memory_order_relaxed) != (size_t)t) {
            size_t workerLastSeq = mNextSeq[t].load(std::memory_order_relaxed) - W;
            if (!anyProcessed || workerLastSeq > lastSeq) {
                lastSeq = workerLastSeq;
                anyProcessed = true;
            }
        }
    }
    size_t offset = anyProcessed ?
        mOffsetRing[lastSeq & (OFFSET_RING_SIZE - 1)].cumulative_offset.load(std::memory_order_relaxed) : 0;
    ftruncate(mFd, offset);
}

void WriterThread::output(){
    if (mPwriteMode) return;  // no-op
    SingleProducerSingleConsumerList<string*>* list = mBufferLists[mWorkingBufferList];
    if(!list->canBeConsumed()) {
        if(mInputCompleted) return;
        // Current slot has no data yet. Block until a producer or
        // setInputCompleted() wakes us via mWriterNotify.
        uint32_t cur = mWriterNotify.load(std::memory_order_acquire);
        hwy::BlockUntilDifferent(cur, mWriterNotify);
        // After wake, return to writerTask loop which re-checks isCompleted()
        return;
    }
    string* str = list->consume();
    mWriter1->write(str->data(), str->length());
    delete str;
    mBufferLength.fetch_sub(1, std::memory_order_release);
    mWorkingBufferList = (mWorkingBufferList+1)%mOptions->thread;
}

void WriterThread::input(int tid, string* data) {
    if (mPwriteMode) {
        inputPwrite(tid, data);
        return;
    }
    mBufferLists[tid]->produce(data);
    mBufferLength.fetch_add(1, std::memory_order_release);
    mWriterNotify.fetch_add(1, std::memory_order_release);
    hwy::WakeAll(mWriterNotify);
}

void WriterThread::inputPwrite(int tid, string* data) {
    size_t bound = libdeflate_gzip_compress_bound(mCompressors[tid], data->size());
    // Grow per-worker buffer if needed
    if (bound > mCompBufSizes[tid]) {
        delete[] mCompBufs[tid];
        mCompBufs[tid] = new char[bound];
        mCompBufSizes[tid] = bound;
    }
    size_t outsize = libdeflate_gzip_compress(mCompressors[tid], data->data(), data->size(),
                                               mCompBufs[tid], bound);
    if (outsize == 0)
        error_exit("libdeflate gzip compression failed");
    delete data;
    const char* writeData = mCompBufs[tid];
    size_t wsize = outsize;

    size_t seq = mNextSeq[tid].load(std::memory_order_relaxed);

    // Wait for previous batch's cumulative offset.
    // Sleep yields CPU to prevent livelock under contention.
    size_t offset = 0;
    if (seq > 0) {
        size_t prevSlot = (seq - 1) & (OFFSET_RING_SIZE - 1);
        while (mOffsetRing[prevSlot].published_seq.load(std::memory_order_acquire) != seq - 1) {
            std::this_thread::sleep_for(std::chrono::microseconds(1));
        }
        offset = mOffsetRing[prevSlot].cumulative_offset.load(std::memory_order_relaxed);
    }

    // Publish offset BEFORE pwrite — next worker starts immediately
    size_t mySlot = seq & (OFFSET_RING_SIZE - 1);
    mOffsetRing[mySlot].cumulative_offset.store(offset + wsize, std::memory_order_relaxed);
    mOffsetRing[mySlot].published_seq.store(seq, std::memory_order_release);

    // pwrite (concurrent with other workers on non-overlapping regions)
    if (wsize > 0) {
        size_t written = 0;
        while (written < wsize) {
            ssize_t ret = pwrite(mFd, writeData + written, wsize - written, offset + written);
            if (ret < 0) {
                if (errno == EINTR) continue;
                error_exit("pwrite failed: " + string(strerror(errno)));
            }
            if (ret == 0)
                error_exit("pwrite returned 0 (disk full?)");
            written += ret;
        }
    }

    // Release store: ensures the pwrite and cumulative_offset publication
    // happen-before the acquire fence in setInputCompletedPwrite().
    mNextSeq[tid].store(
        mNextSeq[tid].load(std::memory_order_relaxed) + mOptions->thread,
        std::memory_order_release);
}

void WriterThread::cleanup() {
    if (mPwriteMode) {
        if (mFd >= 0) { close(mFd); mFd = -1; }
        delete[] mOffsetRing; mOffsetRing = NULL;
        delete[] mNextSeq; mNextSeq = NULL;
        if (mCompressors) {
            for (int t = 0; t < mOptions->thread; t++)
                libdeflate_free_compressor(mCompressors[t]);
            delete[] mCompressors; mCompressors = NULL;
        }
        if (mCompBufs) {
            for (int t = 0; t < mOptions->thread; t++)
                delete[] mCompBufs[t];
            delete[] mCompBufs; mCompBufs = NULL;
        }
        delete[] mCompBufSizes; mCompBufSizes = NULL;
        return;
    }
    deleteWriter();
    if (mBufferLists) {
        for(int t=0; t<mOptions->thread; t++)
            delete mBufferLists[t];
        delete[] mBufferLists;
        mBufferLists = NULL;
    }
}

void WriterThread::deleteWriter() {
    if(mWriter1 != NULL) {
        delete mWriter1;
        mWriter1 = NULL;
    }
}

void WriterThread::initWriter(string filename1, bool isSTDOUT) {
    deleteWriter();
    mWriter1 = new Writer(mOptions, filename1, mOptions->compression, isSTDOUT);
}

void WriterThread::initBufferLists() {
    mBufferLists = new SingleProducerSingleConsumerList<string*>*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++) {
        mBufferLists[t] = new SingleProducerSingleConsumerList<string*>();
    }
}
