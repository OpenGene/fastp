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
    mNextOutputSeq = 0;
    mOrderedOutput = true;

    if (mPwriteMode) {
        mFd = open(mFilename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (mFd < 0)
            error_exit("Failed to open for pwrite: " + mFilename);
        mOffsetRing = new OffsetSlot[OFFSET_RING_SIZE];
        mNextSeq = new size_t[mOptions->thread];
        for (int t = 0; t < mOptions->thread; t++)
            mNextSeq[t] = t;
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
    } else {
        initWriter(filename, isSTDOUT);
        initBufferLists();
        mWorkingBufferList = 0;
        mBufferLength = 0;
    }
}

WriterThread::~WriterThread() {
    cleanup();
}

bool WriterThread::isCompleted()
{
    if (mPwriteMode) return true;  // no writer thread needed
    return mInputCompleted.load(std::memory_order_acquire) && allBufferListsDrained();
}

bool WriterThread::allBufferListsDrained() {
    if(mBufferLists == NULL)
        return true;
    for(int t=0; t<mOptions->thread; t++) {
        if(!mBufferLists[t]->isProducerFinished() || !mBufferLists[t]->isEmpty())
            return false;
    }
    return true;
}

bool WriterThread::setInputCompleted() {
    if (mPwriteMode) {
        setInputCompletedPwrite();
        mInputCompleted.store(true, std::memory_order_release);
        return true;
    }
    for(int t=0; t<mOptions->thread; t++) {
        mBufferLists[t]->setProducerFinished();
    }
    mInputCompleted.store(true, std::memory_order_release);
    return true;
}

void WriterThread::setInputCompletedPwrite() {
    int W = mOptions->thread;
    size_t lastSeq = 0;
    bool anyProcessed = false;
    for (int t = 0; t < W; t++) {
        if (mNextSeq[t] != (size_t)t) {
            size_t workerLastSeq = mNextSeq[t] - W;
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
    if(mOrderedOutput)
        outputOrdered();
    else
        outputReady();
}

void WriterThread::outputOrdered() {
    int tid = mNextOutputSeq % mOptions->thread;
    SingleProducerSingleConsumerList<string*>* list = mBufferLists[tid];
    if(list->canBeConsumed()) {
        string* str = list->consume();
        mWriter1->write(str->data(), str->length());
        delete str;
        mBufferLength--;
        mNextOutputSeq++;
        return;
    }
    usleep(100);
}

void WriterThread::outputReady() {
    for(int i=0; i<mOptions->thread; i++) {
        SingleProducerSingleConsumerList<string*>* list = mBufferLists[mWorkingBufferList];
        if(list->canBeConsumed()) {
            string* str = list->consume();
            mWriter1->write(str->data(), str->length());
            delete str;
            mBufferLength--;
            mWorkingBufferList = (mWorkingBufferList+1)%mOptions->thread;
            return;
        }
        mWorkingBufferList = (mWorkingBufferList+1)%mOptions->thread;
    }
    usleep(100);
}

void WriterThread::input(int tid, string* data) {
    if (mPwriteMode) {
        inputPwrite(tid, data);
        return;
    }
    mBufferLists[tid]->produce(data);
    mBufferLength++;
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

    size_t seq = mNextSeq[tid];

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

    mNextSeq[tid] += mOptions->thread;
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
