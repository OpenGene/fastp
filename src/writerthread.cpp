#include "writerthread.h"
#include "util.h"
#include "common.h"
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
    mMaxPublishedSeq = SIZE_MAX;
    mCompressors = NULL;
    mCompBufs = NULL;
    mCompBufSizes = NULL;
    mOutputRing = NULL;
    mNextExpectedSeq = 0;
    mBufferLength = 0;

    if (mPwriteMode) {
        mFd = open(mFilename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (mFd < 0)
            error_exit("Failed to open for pwrite: " + mFilename);
        mOffsetRing = new OffsetSlot[OFFSET_RING_SIZE];
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
    } else {
        initWriter(filename, isSTDOUT);
        initOutputRing();
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
    mOutputCV.notify_all();
    return true;
}

void WriterThread::setInputCompletedPwrite() {
    size_t maxSeq = mMaxPublishedSeq.load(std::memory_order_acquire);
    size_t offset = (maxSeq == SIZE_MAX) ? 0 :
        mOffsetRing[maxSeq % OFFSET_RING_SIZE].cumulative_offset.load(std::memory_order_relaxed);
    ftruncate(mFd, offset);
}

void WriterThread::output(){
    if (mPwriteMode) return;  // no-op
    size_t want = mNextExpectedSeq.load(std::memory_order_relaxed);
    size_t slot = want % mOutputRingSize;
    if (mOutputRing[slot].seq.load(std::memory_order_acquire) != want) {
        // Wait for input()/setInputCompleted() to notify, with a short timeout
        // as a safety net rather than blind-sleeping every empty check.
        std::unique_lock<std::mutex> lk(mOutputMtx);
        mOutputCV.wait_for(lk, std::chrono::microseconds(100));
        return;
    }
    string* str = mOutputRing[slot].data.load(std::memory_order_relaxed);
    mWriter1->write(str->data(), str->length());
    delete str;
    // Free the slot before advancing so a producer wrapping the ring can't
    // publish into it while it still looks "ready" under the old seq.
    mOutputRing[slot].seq.store(SIZE_MAX, std::memory_order_relaxed);
    mNextExpectedSeq.fetch_add(1, std::memory_order_relaxed);
    mBufferLength--;
}

void WriterThread::input(int tid, size_t seq, string* data) {
    if (mPwriteMode) {
        inputPwrite(tid, seq, data);
        return;
    }
    size_t slot = seq % mOutputRingSize;
    mOutputRing[slot].data.store(data, std::memory_order_relaxed);
    mOutputRing[slot].seq.store(seq, std::memory_order_release);
    mBufferLength++;
    mOutputCV.notify_one();
}

void WriterThread::inputPwrite(int tid, size_t seq, string* data) {
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

    // Wait for previous sequence's cumulative offset.
    // Sleep yields CPU to prevent livelock under contention.
    size_t offset = 0;
    if (seq > 0) {
        size_t prevSlot = (seq - 1) % OFFSET_RING_SIZE;
        while (mOffsetRing[prevSlot].published_seq.load(std::memory_order_acquire) != seq - 1) {
            std::this_thread::sleep_for(std::chrono::microseconds(1));
        }
        offset = mOffsetRing[prevSlot].cumulative_offset.load(std::memory_order_relaxed);
    }

    // Publish offset BEFORE pwrite — whoever holds seq+1 starts immediately
    size_t mySlot = seq % OFFSET_RING_SIZE;
    mOffsetRing[mySlot].cumulative_offset.store(offset + wsize, std::memory_order_relaxed);
    mOffsetRing[mySlot].published_seq.store(seq, std::memory_order_release);

    // Track the highest sequence any worker has published, for truncating
    // the file to the right final size once input is complete (workers no
    // longer publish sequences in a fixed per-worker residue class, so this
    // can't be inferred from any one worker's own progress anymore).
    size_t prevMax = mMaxPublishedSeq.load(std::memory_order_relaxed);
    while ((prevMax == SIZE_MAX || seq > prevMax) &&
           !mMaxPublishedSeq.compare_exchange_weak(prevMax, seq, std::memory_order_relaxed)) {
    }

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
}

void WriterThread::cleanup() {
    if (mPwriteMode) {
        if (mFd >= 0) { close(mFd); mFd = -1; }
        delete[] mOffsetRing; mOffsetRing = NULL;
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
    if (mOutputRing) {
        delete[] mOutputRing;
        mOutputRing = NULL;
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

void WriterThread::initOutputRing() {
    // Must comfortably exceed the largest possible gap between the fastest
    // and slowest worker's progress, which packInMemLimit bounds (mirrors
    // PairEndProcessor::mQueueAssignRingSize's sizing rationale).
    mOutputRingSize = (size_t)(packInMemLimit(mOptions->thread) * 4);
    mOutputRing = new OutputSlot[mOutputRingSize];
}
