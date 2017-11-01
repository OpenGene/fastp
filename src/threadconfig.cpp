#include "threadconfig.h"

ThreadConfig::ThreadConfig(int seqCycles, bool paired){
    mStats1 = new Stats(seqCycles);
    if(paired)
        mStats2 = new Stats(seqCycles);
    else
        mStats2 = NULL;
}