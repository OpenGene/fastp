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

int Stats::statRead(Read* r, char qualifiedQual) {
    int qualified = 0;
    int len = r->length();
    const char* seqstr = r->mSeq.mStr.c_str();
    const char* qualstr = r->mQuality.c_str();

    for(int i=0; i<len; i++) {
        char base = seqstr[i];
        char qual = qualstr[i];
        // get last 3 bits
        char b = base & 0x07;

        const char q20 = '5';
        const char q30 = '?';

        if(qual >= q30) {
            mQ30Bases[b][i]++;
            mQ20Bases[b][i]++;
        } else if(qual >= q20) {
            mQ20Bases[b][i]++;
        }

        mCycleBaseContents[b][i]++;

        if(qual >= qualifiedQual)
            qualified ++;

        mCycleTotalBase[i]++;
        mCycleTotalQual[i]+=qual;

    }

    mReads++;
    return qualified;
}

int Stats::getCycles() {
    return mCycles;
}

void Stats::print() {
    cout << "total reads: " << mReads << endl;
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

