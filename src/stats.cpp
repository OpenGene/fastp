#include "stats.h"

Stats::Stats(int guessedCycles, int bufferMargin){
    mReads = 0;
    mCycles = guessedCycles;

    // extend the buffer to make sure it's long enough
    const int bufLen = guessedCycles + bufferMargin;

    for(int i=0; i<8; i++){
        mQ30Bases[i] = new long[bufLen];
        memset(mQ30Bases[i], 0, sizeof(long) * bufLen);

        mQ20Bases[i] = new long[bufLen];
        memset(mQ20Bases[i], 0, sizeof(long) * bufLen);

        mCycleBaseContents[i] = new long[bufLen];
        memset(mCycleBaseContents[i], 0, sizeof(long) * bufLen);
    }
    mCycleTotalBase = new long[bufLen];
    memset(mCycleTotalBase, 0, sizeof(long)*bufLen);

    mCycleTotalQual = new long[bufLen];
    memset(mCycleTotalQual, 0, sizeof(long)*bufLen);
}

Stats::~Stats() {
    for(int i=0; i<8; i++){
        delete mQ30Bases[i];
        mQ30Bases[i] = NULL;

        delete mQ20Bases[i];
        mQ20Bases[i] = NULL;

        delete mCycleBaseContents[i];
        mCycleBaseContents[i] = NULL;
    }

    delete mCycleTotalBase;
    delete mCycleTotalQual;
}

int Stats::getCycles() {
    return mCycles;
}

Stats* Stats::merge(vector<Stats*>& list) {
    //get the most long cycles
    int cycles = 0;
    for(int t=0; t<list.size(); t++) {
        cycles = max(cycles, list[t]->getCycles());
    }

    Stats* s = new Stats(cycles, 0);

    for(int t=0; t<list.size(); t++) {
        // merge read number
        s->mReads += list[t]->mReads;

        // merge per cycle counting for different bases
        for(int i=0; i<8; i++){
            for(int j=0; j<cycles; j++) {
                s->mQ30Bases[i][j] += list[t]->mQ30Bases[i][j];
                s->mQ20Bases[i][j] += list[t]->mQ20Bases[i][j];
                s->mCycleBaseContents[i][j] += list[t]->mCycleBaseContents[i][j];
            }
        }

        // merge per cycle counting for all bases
        for(int j=0; j<cycles; j++) {
            s->mCycleTotalBase[j] += list[t]->mCycleTotalBase[j];
            s->mCycleTotalQual[j] += list[t]->mCycleTotalQual[j];
        }
    }

    return s;

}

