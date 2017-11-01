#include "threadconfig.h"

ThreadConfig::ThreadConfig(Options* opt, int seqCycles, bool paired){
    mOptions = opt;
    mPreStats1 = new Stats(seqCycles);
    mPreStats2 = new Stats(seqCycles);
    if(paired){
        mPreStats2 = new Stats(seqCycles);
        mPostStats2 = new Stats(seqCycles);
    }
    else {
        mPreStats2 = NULL;
        mPostStats2 = NULL;
    }
}

void ThreadConfig::initWriter(string filename1, string filename2) {
    mWriter1 = new FastqWriter(filename1, mOptions->compression);
    if(filename2 != "" ) {
        mWriter2 = new FastqWriter(filename2, mOptions->compression);
    } else {
        mWriter2 = NULL;
    }
}