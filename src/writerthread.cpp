#include "writerthread.h"
#include "util.h"
#include <memory.h>
#include <unistd.h>

WriterThread::WriterThread(Options* opt, string filename, bool isSTDOUT){
    mOptions = opt;

    mWriter1 = NULL;

    mInputCompleted = false;
    mFilename = filename;

    initWriter(filename, isSTDOUT);
    initBufferLists();
    mWorkingBufferList = 0; // 0 ~ mOptions->thread-1
    mBufferLength = 0;
}

WriterThread::~WriterThread() {
    cleanup();
}

bool WriterThread::isCompleted() 
{
    return mInputCompleted && (mBufferLength==0);
}

bool WriterThread::setInputCompleted() {
    mInputCompleted = true;
    for(int t=0; t<mOptions->thread; t++) {
        mBufferLists[t]->setProducerFinished();
    }
    return true;
}

void WriterThread::output(){
    SingleProducerSingleConsumerList<string*>* list =  mBufferLists[mWorkingBufferList];
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
    mBufferLists[tid]->produce(data);
    mBufferLength++;
}

void WriterThread::cleanup() {
    deleteWriter();
    for(int t=0; t<mOptions->thread; t++) {
        delete mBufferLists[t];
    }
    delete[] mBufferLists;
    mBufferLists = NULL;
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
