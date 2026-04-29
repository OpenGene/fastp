/*
MIT License

Copyright (c) 2021 Shifu Chen <chen@haplox.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

// A ultra-fast lock-free linked list for single-producer, single-consumer threading
// Memory usage overhead: 3M bytes per list, if you want to save memory, please use smaller block and smaller ring buffer
// The type T is usually a pointer, a internal type (such as int, long), or a class supports assignment T a = b;

/* WARNING: only supports up to 1G unconsumed elements in list, 
            which means: produced - consumed must < 1G,
            this is usually much more than enough,
            but if you want to support even more unconsumed elements, 
            please modify the value of blocksRingBufferSize as you want.
*/

#ifndef SINGLE_PRODUCER_SINGLE_CONSUMER_LIST_H
#define SINGLE_PRODUCER_SINGLE_CONSUMER_LIST_H

#include <atomic>
#include <stdio.h>
#include <memory.h>
#include <cassert>
#include <vector>

template<typename T>
struct LockFreeListItem {
public:
    inline LockFreeListItem(T val) {
        value = val;
        nextItemReady = false;
        nextItem = NULL;
    }
    inline LockFreeListItem() {
        nextItem = NULL;
        nextItemReady = false;
    }
    T value;
    LockFreeListItem<T>* nextItem;
    std::atomic_bool nextItemReady;
};

template<typename T>
class SingleProducerSingleConsumerList {
public:
    inline SingleProducerSingleConsumerList() {
        head.store(NULL, std::memory_order_relaxed);
        tail = NULL;
        producerFinished = false;
        consumerFinished = false;
        produced = 0;
        consumed = 0;
        recycled = 0;
        blocksRingBufferSize = 0x1L << 18;
        blocksRingBufferSizeMask = blocksRingBufferSize - 1;
        blocksNum = 0;
        // 2M memory
        blocks = new LockFreeListItem<T>*[blocksRingBufferSize];
        memset(blocks, 0, sizeof(LockFreeListItem<T>*) * blocksRingBufferSize);
    }
    inline ~SingleProducerSingleConsumerList() {
        while(recycled < blocksNum) {
            delete[] blocks[recycled & blocksRingBufferSizeMask];
            blocks[recycled & blocksRingBufferSizeMask] = NULL;
            recycled++;
        }
        delete[] blocks;
        blocks = NULL;
    }
    inline size_t size() {
        return produced.load(std::memory_order_relaxed) - consumed.load(std::memory_order_relaxed);
    }
    inline bool isEmpty() const {
        return head == NULL;
    }
    inline bool canBeConsumed() {
        LockFreeListItem<T>* h = head.load(std::memory_order_acquire);
        if(h == NULL)
            return false;
        // `nextItemReady` is a publication barrier for `nextItem`.
        // The last node has no successor, so `nextItemReady` may remain false;
        // it must still be consumable to avoid writer stalls when many queues exist.
        return h->nextItemReady.load(std::memory_order_acquire)
            || (h == tail)
            || producerFinished.load(std::memory_order_acquire);
    }
    inline void produce(T val) {
        LockFreeListItem<T>* item = makeItem(val);
        if(head.load(std::memory_order_relaxed) == NULL) {
            tail = item;
            // Release store: publishing head to consumer thread.
            // All writes to *item are ordered before this store.
            head.store(item, std::memory_order_release);
            // Signal the first item is consumable (no predecessor to set this)
            item->nextItemReady.store(true, std::memory_order_release);
        } else {
            tail->nextItem = item;
            // Release store: ensures nextItem write visible before nextItemReady flag.
            tail->nextItemReady.store(true, std::memory_order_release);
            tail = item;
        }
        produced.fetch_add(1, std::memory_order_relaxed);
    }
    inline T consume() {
        LockFreeListItem<T>* h = head.load(std::memory_order_acquire);
        assert(h != NULL);
        T val = h->value;
        // Advance head; release so next canBeConsumed() acquire sees updated state.
        head.store(h->nextItem, std::memory_order_release);
        unsigned long _c = consumed.fetch_add(1, std::memory_order_relaxed) + 1;
        if((_c & 0xFFF) == 0)
            recycle();
        return val;
    }
    inline bool isProducerFinished() {
        return producerFinished.load(std::memory_order_acquire);
    }
    inline bool isConsumerFinished() {
        return consumerFinished.load(std::memory_order_acquire);
    }
    inline void setProducerFinished() {
        producerFinished.store(true, std::memory_order_release);
    }
    inline void setConsumerFinished() {
        consumerFinished.store(true, std::memory_order_release);
    }
private:
    // blockized list
    inline LockFreeListItem<T>* makeItem(T val) {
        unsigned long _p = produced.load(std::memory_order_relaxed);
        unsigned long blk = _p >> 12;
        unsigned long idx = _p & 0xFFF;
        size_t size = 0x01<<12;
        if(blocksNum <= blk) {
            LockFreeListItem<T>* buffer = new LockFreeListItem<T>[size];
            memset(buffer, 0, sizeof(LockFreeListItem<T>) * size);
            blocks[blocksNum & blocksRingBufferSizeMask] = buffer;
            blocksNum++;
        }
        LockFreeListItem<T>* item = blocks[blk & blocksRingBufferSizeMask]+idx;
        item->value = val;
        return item;
    }

    inline void recycle() {
        unsigned long blk = consumed.load(std::memory_order_relaxed) >> 12;
        while((recycled+1) < blk) {
            delete[] blocks[recycled & blocksRingBufferSizeMask];
            blocks[recycled & blocksRingBufferSizeMask] = NULL;
            recycled++;
        }
    }

private:
    std::atomic<LockFreeListItem<T>*> head;
    LockFreeListItem<T>* tail;       // tail is producer-private, no atomic needed
    LockFreeListItem<T>** blocks;
    std::atomic_bool producerFinished;
    std::atomic_bool consumerFinished;
    std::atomic<unsigned long> produced;
    std::atomic<unsigned long> consumed;
    unsigned long recycled;
    unsigned long blocksRingBufferSize;
    unsigned long blocksRingBufferSizeMask;
    unsigned long blocksNum;
};

#endif
