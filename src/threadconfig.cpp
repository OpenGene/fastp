#include "threadconfig.h"

ThreadConfig::ThreadConfig(Options* opt, int seqCycles, bool paired){
    mOptions = opt;
    mPreStats1 = new Stats(seqCycles);
    mPostStats1 = new Stats(seqCycles);
    if(paired){
        mPreStats2 = new Stats(seqCycles);
        mPostStats2 = new Stats(seqCycles);
    }
    else {
        mPreStats2 = NULL;
        mPostStats2 = NULL;
    }
    mWriter1 = NULL;
    mWriter2 = NULL;
}

ThreadConfig::~ThreadConfig() {
    if(mWriter1 != NULL) {
        delete mWriter1;
        mWriter1 = NULL;
    }
    if(mWriter2 != NULL) {
        delete mWriter2;
        mWriter2 = NULL;
    }
}

void ThreadConfig::initWriter(string filename1) {
    mWriter1 = new Writer(filename1, mOptions->compression);
}

void ThreadConfig::initWriter(string filename1, string filename2) {
    mWriter1 = new Writer(filename1, mOptions->compression);
    mWriter2 = new Writer(filename2, mOptions->compression);
}

void ThreadConfig::initWriter(ofstream* stream) {
    mWriter1 = new Writer(stream);
}

void ThreadConfig::initWriter(ofstream* stream1, ofstream* stream2) {

    mWriter1 = new Writer(stream1);
    mWriter2 = new Writer(stream2);
}

void ThreadConfig::initWriter(gzFile gzfile) {

    mWriter1 = new Writer(gzfile);
}

void ThreadConfig::initWriter(gzFile gzfile1, gzFile gzfile2) {

    mWriter1 = new Writer(gzfile1);
    mWriter2 = new Writer(gzfile2);
}