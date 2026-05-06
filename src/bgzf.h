#ifndef FASTP_BGZF_H
#define FASTP_BGZF_H

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <atomic>
#include <thread>
#include <vector>
#include <isa-l/igzip_lib.h>
#include <hwy/contrib/thread_pool/futex.h>

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
// Slot lifecycle: FREE -> COMPRESSED -> DECOMPRESSING -> READY -> FREE
// Synchronization: Highway futex (BlockUntilDifferent/WakeAll) on per-slot
// atomic state, plus a global mSlotNotify counter for decompressor wakeup.
class BgzfMtReader {
    static const uint32_t STATE_FREE          = 0;
    static const uint32_t STATE_COMPRESSED    = 1;
    static const uint32_t STATE_DECOMPRESSING = 2;
    static const uint32_t STATE_READY         = 3;
    static const uint32_t STATE_DONE          = 4;

    struct alignas(64) Slot {
        unsigned char comp[BGZF_MAX_BLOCK_SIZE];
        unsigned char decomp[BGZF_MAX_BLOCK_SIZE];
        int compLen;
        int decompLen;
        std::atomic<uint32_t> state{STATE_FREE};
    };

public:
    BgzfMtReader(FILE* fp, int threadBudget = 0)
        : mFp(fp), mConsumeIdx(0), mConsumeOffset(0),
          mProduceIdx(0), mStop(false) {
        mSlotNotify.store(0, std::memory_order_relaxed);
        int cpus = std::thread::hardware_concurrency();
        if (cpus < 2) cpus = 2;
        mMaxPool = (threadBudget > 0) ? threadBudget : std::max(1, cpus / 2);
        mRingSize = mMaxPool * 4;
        if (mRingSize < 16) mRingSize = 16;
        if (mRingSize > 64) mRingSize = 64;
        mSlots = new Slot[mRingSize];

        for (int i = 0; i < mMaxPool; i++)
            mPool.push_back(new std::thread(&BgzfMtReader::decompWorker, this));
        mReaderThread = new std::thread(&BgzfMtReader::readerLoop, this);
    }

    ~BgzfMtReader() {
        mStop.store(true, std::memory_order_release);
        // Wake all waiters so they see mStop
        notifyAll();
        if (mReaderThread) { mReaderThread->join(); delete mReaderThread; }
        for (auto* t : mPool) { t->join(); delete t; }
        delete[] mSlots;
    }

    int read(char* outBuf, int outBufSize) {
        int filled = 0;
        while (filled < outBufSize) {
            Slot& s = mSlots[mConsumeIdx % mRingSize];
            // Wait for slot to become READY or DONE
            while (true) {
                uint32_t v = s.state.load(std::memory_order_acquire);
                if (v == STATE_READY || v == STATE_DONE) break;
                hwy::BlockUntilDifferent(v, s.state);
            }
            if (s.state.load(std::memory_order_acquire) == STATE_DONE) break;

            int avail = s.decompLen - mConsumeOffset;
            int tocopy = (outBufSize - filled) < avail ? (outBufSize - filled) : avail;
            memcpy(outBuf + filled, s.decomp + mConsumeOffset, tocopy);
            filled += tocopy;
            mConsumeOffset += tocopy;

            if (mConsumeOffset >= s.decompLen) {
                s.state.store(STATE_FREE, std::memory_order_release);
                hwy::WakeAll(s.state);
                mConsumeOffset = 0;
                mConsumeIdx++;
                // Notify producer that a slot is free
                mSlotNotify.fetch_add(1, std::memory_order_release);
                hwy::WakeAll(mSlotNotify);
            }
        }
        return filled;
    }

private:
    void notifyAll() {
        // Wake all threads blocked on slot states or mSlotNotify
        for (int i = 0; i < mRingSize; i++) {
            hwy::WakeAll(mSlots[i].state);
        }
        mSlotNotify.fetch_add(1, std::memory_order_release);
        hwy::WakeAll(mSlotNotify);
    }

    void readerLoop() {
        unsigned char header[BGZF_HEADER_SIZE];

        while (!mStop.load(std::memory_order_acquire)) {
            Slot& s = mSlots[mProduceIdx % mRingSize];
            // Wait for slot to become FREE
            while (true) {
                if (mStop.load(std::memory_order_acquire)) return;
                uint32_t v = s.state.load(std::memory_order_acquire);
                if (v == STATE_FREE) break;
                hwy::BlockUntilDifferent(v, s.state);
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
            s.state.store(STATE_COMPRESSED, std::memory_order_release);
            hwy::WakeAll(s.state);
            mProduceIdx++;
            // Notify decompressors that a slot has data
            mSlotNotify.fetch_add(1, std::memory_order_release);
            hwy::WakeAll(mSlotNotify);
        }
    }

    void decompWorker() {
        while (!mStop.load(std::memory_order_acquire)) {
            Slot* target = nullptr;
            if (!claimSlot(&target)) {
                if (mStop.load(std::memory_order_acquire)) return;
                // No slot available — block until state changes somewhere
                uint32_t cur = mSlotNotify.load(std::memory_order_acquire);
                hwy::BlockUntilDifferent(cur, mSlotNotify);
                continue;
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

            target->state.store(STATE_READY, std::memory_order_release);
            hwy::WakeAll(target->state);
            // Notify consumer that data is ready
            mSlotNotify.fetch_add(1, std::memory_order_release);
            hwy::WakeAll(mSlotNotify);
        }
    }

    bool claimSlot(Slot** out) {
        int head = mConsumeIdx;
        int tail = mProduceIdx;
        for (int i = head; i < tail; i++) {
            Slot& s = mSlots[i % mRingSize];
            uint32_t expected = STATE_COMPRESSED;
            if (s.state.compare_exchange_strong(expected, STATE_DECOMPRESSING,
                    std::memory_order_acq_rel)) {
                *out = &s;
                return true;
            }
        }
        return false;
    }

    void markDone(Slot& s) {
        s.state.store(STATE_DONE, std::memory_order_release);
        hwy::WakeAll(s.state);
        mSlotNotify.fetch_add(1, std::memory_order_release);
        hwy::WakeAll(mSlotNotify);
    }

    FILE* mFp;
    Slot* mSlots;
    int mRingSize;
    std::atomic<int> mConsumeIdx;
    int mConsumeOffset;
    std::atomic<int> mProduceIdx;
    int mMaxPool;
    std::atomic<bool> mStop;
    std::atomic<uint32_t> mSlotNotify;  // monotonic counter for cross-thread wakeup

    std::thread* mReaderThread;
    std::vector<std::thread*> mPool;
};

#endif
