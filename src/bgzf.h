#ifndef FASTP_BGZF_H
#define FASTP_BGZF_H

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <isa-l/igzip_lib.h>

static const int BGZF_HEADER_SIZE = 18;
static const int BGZF_MAX_BLOCK_SIZE = 65536;

static inline bool isBgzf(FILE* fp) {
    unsigned char h[18];
    long pos = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    size_t n = fread(h, 1, 18, fp);
    fseek(fp, pos, SEEK_SET);
    if (n < 18) return false;
    return h[0]==0x1f && h[1]==0x8b && h[2]==0x08 && (h[3]&0x04)
        && h[12]==0x42 && h[13]==0x43 && (h[14]|(h[15]<<8))==2;
}

static inline uint32_t bgzfBlockSize(const unsigned char* h) {
    if (h[0]!=0x1f || h[1]!=0x8b) return 0;
    return (uint32_t)(h[16] | (h[17]<<8)) + 1;
}

// Parallel BGZF decompression with fixed thread pool.
// All decompress threads start upfront; ring buffer provides backpressure.
// Slot lifecycle: FREE → COMPRESSED → DECOMPRESSING → READY → FREE
class BgzfMtReader {
    enum SlotState { FREE, COMPRESSED, DECOMPRESSING, READY, DONE };

    struct alignas(64) Slot {
        unsigned char comp[BGZF_MAX_BLOCK_SIZE];
        unsigned char decomp[BGZF_MAX_BLOCK_SIZE];
        int compLen;
        int decompLen;
        std::atomic<SlotState> state{FREE};
    };

public:
    // threadBudget: max decompress threads allowed for this reader.
    // 0 = auto (use half of available CPU cores).
    // Caller should compute: (hardware_concurrency - workers - readers - writers) / num_gz_inputs
    BgzfMtReader(FILE* fp, int threadBudget = 0)
        : mFp(fp), mConsumeIdx(0), mConsumeOffset(0),
          mProduceIdx(0), mStop(false) {
        int cpus = std::thread::hardware_concurrency();
        if (cpus < 2) cpus = 2;
        if (cpus > 8) cpus = 8; // cap at 8 cores for now, to avoid starving workers
        mMaxPool = (threadBudget > 0) ? threadBudget : std::max(1, cpus / 2);
        mRingSize = mMaxPool * 4;
        if (mRingSize < 16) mRingSize = 16;
        if (mRingSize > 64) mRingSize = 64;  // cap at ~8MB per reader
        mSlots = new Slot[mRingSize];

        // Start all decompress threads upfront — ring buffer handles backpressure
        for (int i = 0; i < mMaxPool; i++)
            mPool.push_back(new std::thread(&BgzfMtReader::decompWorker, this));
        mReaderThread = new std::thread(&BgzfMtReader::readerLoop, this);
    }

    ~BgzfMtReader() {
        // Publish mStop while holding every CV's mutex before notifying. mStop is a
        // wait predicate for all three CVs; a thread that has evaluated its predicate
        // as false but not yet parked still holds that mutex, so taking all three here
        // serializes the stop+notify against that gap. Storing mStop outside the locks
        // (as before) lets notify_all() slip into that window and be missed, leaving a
        // worker parked forever and the join() below hung — an intermittent deadlock.
        {
            std::lock_guard<std::mutex> l1(mProduceMtx);
            std::lock_guard<std::mutex> l2(mDecompMtx);
            std::lock_guard<std::mutex> l3(mConsumeMtx);
            mStop = true;
        }
        mDecompCv.notify_all();
        mProduceCv.notify_all();
        mConsumeCv.notify_all();
        if (mReaderThread) { mReaderThread->join(); delete mReaderThread; }
        for (auto* t : mPool) { t->join(); delete t; }
        delete[] mSlots;
    }

    int read(char* outBuf, int outBufSize) {
        int filled = 0;
        while (filled < outBufSize) {
            Slot& s = mSlots[mConsumeIdx % mRingSize];
            {
                std::unique_lock<std::mutex> lk(mConsumeMtx);
                mConsumeCv.wait(lk, [&]() {
                    SlotState v = s.state.load(std::memory_order_acquire);
                    return v == READY || v == DONE;
                });
            }
            if (s.state.load(std::memory_order_acquire) == DONE) break;

            int avail = s.decompLen - mConsumeOffset;
            int tocopy = (outBufSize - filled) < avail ? (outBufSize - filled) : avail;
            memcpy(outBuf + filled, s.decomp + mConsumeOffset, tocopy);
            filled += tocopy;
            mConsumeOffset += tocopy;

            if (mConsumeOffset >= s.decompLen) {
                // Publish FREE under mProduceMtx so readerLoop, which waits on
                // mProduceCv for this slot to free up, cannot miss the wakeup.
                {
                    std::lock_guard<std::mutex> lk(mProduceMtx);
                    s.state.store(FREE, std::memory_order_release);
                }
                mConsumeOffset = 0;
                mConsumeIdx++;
                mProduceCv.notify_one();
            }
        }
        return filled;
    }

private:
    void readerLoop() {
        unsigned char header[BGZF_HEADER_SIZE];

        while (!mStop) {
            Slot& s = mSlots[mProduceIdx % mRingSize];
            {
                std::unique_lock<std::mutex> lk(mProduceMtx);
                mProduceCv.wait(lk, [&]() {
                    return s.state.load(std::memory_order_acquire) == FREE || mStop;
                });
                if (mStop) break;
            }

            size_t n = fread(header, 1, BGZF_HEADER_SIZE, mFp);
            if (n < BGZF_HEADER_SIZE) { markDone(s); break; }

            uint32_t bsize = bgzfBlockSize(header);
            if (bsize == 0 || bsize > BGZF_MAX_BLOCK_SIZE) { markDone(s); break; }

            memcpy(s.comp, header, BGZF_HEADER_SIZE);
            size_t rest = bsize - BGZF_HEADER_SIZE;
            if (fread(s.comp + BGZF_HEADER_SIZE, 1, rest, mFp) < rest) { markDone(s); break; }

            if (bsize == 28) {
                uint32_t isize = s.comp[24]|(s.comp[25]<<8)|(s.comp[26]<<16)|(s.comp[27]<<24);
                if (isize == 0) { markDone(s); break; }
            }

            s.compLen = bsize;
            // Publish the COMPRESSED slot (state + mProduceIdx, which claimSlot scans)
            // under mDecompMtx so a decompWorker parked on mDecompCv cannot miss it.
            {
                std::lock_guard<std::mutex> lk(mDecompMtx);
                s.state.store(COMPRESSED, std::memory_order_release);
                mProduceIdx++;
            }
            mDecompCv.notify_all();
        }
    }

    void decompWorker() {
        while (!mStop) {
            Slot* target = nullptr;
            {
                std::unique_lock<std::mutex> lk(mDecompMtx);
                mDecompCv.wait(lk, [&]() {
                    return mStop.load(std::memory_order_relaxed) || claimSlot(&target);
                });
                if (mStop && !target) return;
                if (!target) continue;
            }

            struct inflate_state ist;
            isal_inflate_init(&ist);
            ist.crc_flag = ISAL_GZIP;
            ist.next_in = target->comp;
            ist.avail_in = target->compLen;
            ist.next_out = target->decomp;
            ist.avail_out = BGZF_MAX_BLOCK_SIZE;
            int ret = isal_inflate_stateless(&ist);
            target->decompLen = (ret == ISAL_DECOMP_OK) ? (int)ist.total_out : 0;

            // Publish READY under mConsumeMtx so the consumer parked in read() on
            // mConsumeCv cannot miss the wakeup for this slot.
            {
                std::lock_guard<std::mutex> lk(mConsumeMtx);
                target->state.store(READY, std::memory_order_release);
            }
            mConsumeCv.notify_one();
        }
    }

    bool claimSlot(Slot** out) {
        int head = mConsumeIdx;
        int tail = mProduceIdx;
        for (int i = head; i < tail; i++) {
            Slot& s = mSlots[i % mRingSize];
            SlotState expected = COMPRESSED;
            if (s.state.compare_exchange_strong(expected, DECOMPRESSING,
                    std::memory_order_acq_rel)) {
                *out = &s;
                return true;
            }
        }
        return false;
    }

    void markDone(Slot& s) {
        // Publish DONE under mConsumeMtx before notifying: the consumer in read()
        // waits on mConsumeCv for READY/DONE, and at EOF this is the terminal
        // notification (no further producer follows), so a lost wakeup here would
        // strand the reader thread in read() forever.
        {
            std::lock_guard<std::mutex> lk(mConsumeMtx);
            s.state.store(DONE, std::memory_order_release);
        }
        mConsumeCv.notify_all();
        mDecompCv.notify_all();
    }

    FILE* mFp;  // not owned — caller must keep open for lifetime of BgzfMtReader
    Slot* mSlots;
    int mRingSize;
    std::atomic<int> mConsumeIdx;
    int mConsumeOffset;  // only accessed by consumer thread
    std::atomic<int> mProduceIdx;
    int mMaxPool;
    std::atomic<bool> mStop;

    std::thread* mReaderThread;
    std::vector<std::thread*> mPool;

    std::mutex mConsumeMtx, mProduceMtx, mDecompMtx;
    std::condition_variable mConsumeCv, mProduceCv, mDecompCv;
};

#endif
