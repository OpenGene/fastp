#include "writerthread.h"
#include "util.h"
#include <memory.h>
#include <unistd.h>
#include <fcntl.h>
#include <cerrno>
#include <cstring>
#include <thread>

WriterThread::WriterThread(Options* opt, string filename, bool isSTDOUT){
    mOptions = opt;
    mWriter1 = NULL;
    mInputCompleted = false;
    mFilename = filename;

    mPwriteMode = !isSTDOUT && mOptions->thread > 1;
    mCompressInPwrite = mPwriteMode && ends_with(filename, ".gz");
    mFd = -1;
    mOffsetRing = NULL;
    mNextSeq = NULL;
    mCompressors = NULL;
    mCompBufs = NULL;
    mCompBufSize = 0;
    mOrderedMode = false;
    mOrderedTotal = 0;
    mOrderedRing = NULL;
    mOrderedWriteCursor = 0;
    mBufferLists = NULL;

    if (mPwriteMode) {
        mFd = open(mFilename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (mFd < 0)
            error_exit("Failed to open for pwrite: " + mFilename);
        mOffsetRing = new OffsetSlot[OFFSET_RING_SIZE];
        mNextSeq = new size_t[mOptions->thread];
        for (int t = 0; t < mOptions->thread; t++)
            mNextSeq[t] = t;
        if (mCompressInPwrite) {
            mCompressors = new libdeflate_compressor*[mOptions->thread];
            for (int t = 0; t < mOptions->thread; t++)
                mCompressors[t] = libdeflate_alloc_compressor(mOptions->compression);
            mCompBufSize = PACK_SIZE * 500;  // ~500 bytes/read worst case
            mCompBufs = new char*[mOptions->thread];
            for (int t = 0; t < mOptions->thread; t++)
                mCompBufs[t] = new char[mCompBufSize];
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
    if (mOrderedMode)
        return mInputCompleted && mOrderedWriteCursor.load(std::memory_order_acquire) >= mOrderedTotal;
    return mInputCompleted && (mBufferLength==0);
}

bool WriterThread::setInputCompleted() {
    if (mPwriteMode) {
        if (mOrderedMode) {
            if (mOrderedTotal > 0) {
                size_t lastSeq = mOrderedTotal - 1;
                size_t lastSlot = lastSeq & (OFFSET_RING_SIZE - 1);
                for (int spins = 0; mOffsetRing[lastSlot].published_seq.load(std::memory_order_acquire) != lastSeq; ) {
                    if (++spins > 256) {
                        usleep(1);
                        spins = 0;
                    } else {
#if defined(__aarch64__)
                        __asm__ volatile("yield");
#elif defined(__x86_64__) || defined(__i386__)
                        __asm__ volatile("pause");
#endif
                    }
                }
                size_t offset = mOffsetRing[lastSlot].cumulative_offset.load(std::memory_order_relaxed);
                ftruncate(mFd, offset);
            } else {
                ftruncate(mFd, 0);
            }
        } else {
            setInputCompletedPwrite();
        }
        mInputCompleted = true;
        return true;
    }
    mInputCompleted = true;
    if (!mOrderedMode) {
        for(int t=0; t<mOptions->thread; t++) {
            mBufferLists[t]->setProducerFinished();
        }
    }
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

void WriterThread::setOrderedMode(size_t totalPacks) {
    mOrderedMode = true;
    mOrderedTotal = totalPacks;
    if (!mPwriteMode) {
        mOrderedRing = new std::atomic<string*>[OFFSET_RING_SIZE];
        for (int i = 0; i < OFFSET_RING_SIZE; i++)
            mOrderedRing[i].store(nullptr, std::memory_order_relaxed);
    }
}

void WriterThread::output(){
    if (mPwriteMode) return;  // no-op
    if (mOrderedMode) {
        size_t cursor = mOrderedWriteCursor.load(std::memory_order_relaxed);
        if (cursor >= mOrderedTotal) return;
        size_t slot = cursor & (OFFSET_RING_SIZE - 1);
        string* str = mOrderedRing[slot].load(std::memory_order_acquire);
        if (str) {
            mWriter1->write(str->data(), str->length());
            delete str;
            mOrderedRing[slot].store(nullptr, std::memory_order_release);
            mOrderedWriteCursor.store(cursor + 1, std::memory_order_release);
            mBufferLength--;
        } else {
            usleep(100);
        }
        return;
    }
    SingleProducerSingleConsumerList<string*>* list = mBufferLists[mWorkingBufferList];
    if(!list->canBeConsumed()) {
        usleep(100);
    } else {
        string* str = list->consume();
        mWriter1->write(str->data(), str->length());
        delete str;
        mBufferLength--;
        mWorkingBufferList = (mWorkingBufferList+1)%mOptions->thread;
    }
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
    size_t seq = mNextSeq[tid];
    doInputPwrite(tid, data, seq);
    mNextSeq[tid] += mOptions->thread;
}

void WriterThread::doInputPwrite(int tid, string* data, size_t seq) {
    const char* writeData;
    size_t wsize;

    if (mCompressInPwrite) {
        size_t bound = libdeflate_gzip_compress_bound(mCompressors[tid], data->size());
        if (bound > mCompBufSize) {
            delete[] mCompBufs[tid];
            mCompBufs[tid] = new char[bound];
        }
        size_t outsize = libdeflate_gzip_compress(mCompressors[tid], data->data(), data->size(),
                                                   mCompBufs[tid], bound);
        if (outsize == 0)
            error_exit("libdeflate gzip compression failed");
        delete data;
        writeData = mCompBufs[tid];
        wsize = outsize;
    } else {
        writeData = data->data();
        wsize = data->size();
    }

    // Wait for previous batch's cumulative offset
    size_t offset = 0;
    if (seq > 0) {
        size_t prevSlot = (seq - 1) & (OFFSET_RING_SIZE - 1);
        for (int spins = 0; mOffsetRing[prevSlot].published_seq.load(std::memory_order_acquire) != seq - 1; ) {
            if (++spins > 256) {
                usleep(1);
                spins = 0;
            } else {
#if defined(__aarch64__)
                __asm__ volatile("yield");
#elif defined(__x86_64__) || defined(__i386__)
                __asm__ volatile("pause");
#endif
            }
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

    if (!mCompressInPwrite)
        delete data;
}

void WriterThread::inputOrderedRing(string* data, size_t seq) {
    while (seq - mOrderedWriteCursor.load(std::memory_order_acquire) >= (size_t)OFFSET_RING_SIZE) {
        usleep(10);
    }
    size_t slot = seq & (OFFSET_RING_SIZE - 1);
    mOrderedRing[slot].store(data, std::memory_order_release);
    mBufferLength++;
}

void WriterThread::inputWithSeq(int tid, string* data, size_t seq) {
    if (mPwriteMode) {
        doInputPwrite(tid, data, seq);
    } else {
        inputOrderedRing(data, seq);
    }
}

void WriterThread::cleanup() {
    if (mOrderedRing) {
        for (int i = 0; i < OFFSET_RING_SIZE; i++) {
            string* s = mOrderedRing[i].load(std::memory_order_relaxed);
            if (s) delete s;
        }
        delete[] mOrderedRing;
        mOrderedRing = NULL;
    }
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
