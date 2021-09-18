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

// A lock-free linked list for single-producer, single-consumer threading

#ifndef SINGLEPRODUCERSINGLECONSUMERLIST_H
#define SINGLEPRODUCERSINGLECONSUMERLIST_H

#include <atomic>
#include <stdio.h>
#include <memory.h>
#include <cassert>

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
    }
    inline bool canBeConsumed() {
        if(head == NULL)
            return false;
        return head->nextItemReady || producerFinished;
    }
    inline void produce(T val) {
        LockFreeListItem<T>* item = new LockFreeListItem<T>(val);
        if(head==NULL) {
            head = item;
            tail = item;
        } else {
            tail->nextItem = item;
            tail->nextItemReady = true;
            tail = item;
        }
    }
    inline T consume() {
        assert(head != NULL);
        T val = head->value;
        LockFreeListItem<T>* tmp = head;
        head = head->nextItem;
        delete tmp;
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
    LockFreeListItem<T>* head;
    LockFreeListItem<T>* tail;
    std::atomic_bool producerFinished;
    std::atomic_bool consumerFinished;
};

#endif
