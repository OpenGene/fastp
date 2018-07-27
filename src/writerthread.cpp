#include "writerthread.h"
#include "util.h"
#include <memory.h>
#include <unistd.h>

WriterThread::WriterThread(Options* opt, string filename){
    mOptions = opt;

    mWriter1 = NULL;

    mInputCounter = 0;
    mOutputCounter = 0;
    mInputCompleted = false;
    mFilename = filename;

    mRingBuffer = new char*[PACK_NUM_LIMIT];
    memset(mRingBuffer, 0, sizeof(char*) * PACK_NUM_LIMIT);
    mRingBufferSizes = new size_t[PACK_NUM_LIMIT];
    memset(mRingBufferSizes, 0, sizeof(size_t) * PACK_NUM_LIMIT);
    initWriter(filename);
}

WriterThread::~WriterThread() {
    cleanup();
    delete mRingBuffer;
}

bool WriterThread::isCompleted() 
{
    return mInputCompleted && (mOutputCounter == mInputCounter);
}

bool WriterThread::setInputCompleted() {
    mInputCompleted = true;
    return true;
}

void WriterThread::output(){
    if(mOutputCounter >= mInputCounter) {
        usleep(100);
    }
    while( mOutputCounter < mInputCounter) 
    {
        mWriter1->write(mRingBuffer[mOutputCounter], mRingBufferSizes[mOutputCounter]);
        delete mRingBuffer[mOutputCounter];
        mRingBuffer[mOutputCounter] = NULL;
        mOutputCounter++;
    }
}

void  WriterThread::input(char* data, size_t size){
    mRingBuffer[mInputCounter] = data;
    mRingBufferSizes[mInputCounter] = size;
    mInputCounter++;
}

void WriterThread::cleanup() {
    deleteWriter();
}

void WriterThread::deleteWriter() {
    if(mWriter1 != NULL) {
        delete mWriter1;
        mWriter1 = NULL;
    }
}

void WriterThread::initWriter(string filename1) {
    deleteWriter();
    mWriter1 = new Writer(filename1, mOptions->compression);
}

void WriterThread::initWriter(ofstream* stream) {
    deleteWriter();
    mWriter1 = new Writer(stream);
}

void WriterThread::initWriter(gzFile gzfile) {
    deleteWriter();
    mWriter1 = new Writer(gzfile);
}

long WriterThread::bufferLength(){
    return mInputCounter - mOutputCounter;
}
