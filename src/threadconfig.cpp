#include "threadconfig.h"

ThreadConfig::ThreadConfig(int seqCycles, bool paired){
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