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
        head = NULL;
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
        return produced -  consumed;
    }
    inline bool canBeConsumed() {
        if(head == NULL)
            return false;
        return head->nextItemReady || producerFinished;
    }
    inline void produce(T val) {
        LockFreeListItem<T>* item = makeItem(val);
        if(head==NULL) {
            head = item;
            tail = item;
        } else {
            tail->nextItem = item;
            tail->nextItemReady = true;
            tail = item;
        }
        produced++;
    }
    inline T consume() {
        assert(head != NULL);
        T val = head->value;
        head = head->nextItem;
        consumed++;
        if((consumed & 0xFFF) == 0)
            recycle();
        return val;
    }
    inline bool isProducerFinished() {
        return producerFinished;
    }
    inline bool isConsumerFinished() {
        return consumerFinished;
    }
    inline void setProducerFinished() {
        producerFinished = true;
    }
    inline void setConsumerFinished() {
        consumerFinished = true;
    }
private:
    // blockized list
    inline LockFreeListItem<T>* makeItem(T val) {
        unsigned long blk = produced >> 12;
        unsigned long idx = produced & 0xFFF;
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
        unsigned long blk = consumed >> 12;
        while((recycled+1) < blk) {
            delete[] blocks[recycled & blocksRingBufferSizeMask];
            blocks[recycled & blocksRingBufferSizeMask] = NULL;
            recycled++;
        }
    }

private:
    LockFreeListItem<T>* head;
    LockFreeListItem<T>* tail;
    LockFreeListItem<T>** blocks;
    std::atomic_bool producerFinished;
    std::atomic_bool consumerFinished;
    unsigned long produced;
    unsigned long consumed;
    unsigned long recycled;
    unsigned long blocksRingBufferSize;
    unsigned long blocksRingBufferSizeMask;
    unsigned long blocksNum;
};

#endif
