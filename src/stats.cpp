#include "stats.h"

Stats::Stats(int guessedCycles, int bufferMargin){
    mReads = 0;
    mCycles = guessedCycles;
    mBases = 0;
    mQ20Total = 0;
    mQ30Total = 0;
    summarized = false;

    // extend the buffer to make sure it's long enough
    mBufLen = guessedCycles + bufferMargin;

    for(int i=0; i<8; i++){
        mQ20Bases[i] = 0;
        mQ30Bases[i] = 0;

        mCycleQ30Bases[i] = new long[mBufLen];
        memset(mCycleQ30Bases[i], 0, sizeof(long) * mBufLen);

        mCycleQ20Bases[i] = new long[mBufLen];
        memset(mCycleQ20Bases[i], 0, sizeof(long) * mBufLen);

        mCycleBaseContents[i] = new long[mBufLen];
        memset(mCycleBaseContents[i], 0, sizeof(long) * mBufLen);
    }
    mCycleTotalBase = new long[mBufLen];
    memset(mCycleTotalBase, 0, sizeof(long)*mBufLen);

    mCycleTotalQual = new long[mBufLen];
    memset(mCycleTotalQual, 0, sizeof(long)*mBufLen);
}

Stats::~Stats() {
    for(int i=0; i<8; i++){
        delete mCycleQ30Bases[i];
        mCycleQ30Bases[i] = NULL;

        delete mCycleQ20Bases[i];
        mCycleQ20Bases[i] = NULL;

        delete mCycleBaseContents[i];
        mCycleBaseContents[i] = NULL;
    }

    delete mCycleTotalBase;
    delete mCycleTotalQual;
}

void Stats::summarize() {

    // first get the cycle and count total bases
    for(int c=0; c<mBufLen; c++) {
        mBases += mCycleTotalBase[c];
        if (mCycleTotalBase[c] == 0){
            mCycles = c;
            break;
        }
    }

    // Q20 Q30
    for(int i=0; i<8; i++) {
        for(int c=0; c<mCycles; c++) {
            mQ20Bases[i] += mCycleQ20Bases[i][c];
            mQ30Bases[i] += mCycleQ30Bases[i][c];
        }
        mQ20Total += mQ20Bases[i];
        mQ30Total += mQ30Bases[i];
    }
    summarized = true;
}

void Stats::statRead(Read* r, int& lowQualNum, int& nBaseNum, char qualifiedQual) {
    lowQualNum = 0;
    nBaseNum = 0;
    int len = r->length();
    const char* seqstr = r->mSeq.mStr.c_str();
    const char* qualstr = r->mQuality.c_str();

    for(int i=0; i<len; i++) {
        char base = seqstr[i];
        char qual = qualstr[i];
        // get last 3 bits
        char b = base & 0x07;

        if(base == 'N')
            nBaseNum++;

        const char q20 = '5';
        const char q30 = '?';

        if(qual >= q30) {
            mCycleQ30Bases[b][i]++;
            mCycleQ20Bases[b][i]++;
        } else if(qual >= q20) {
            mCycleQ20Bases[b][i]++;
        }

        mCycleBaseContents[b][i]++;

        if(qual < qualifiedQual)
            lowQualNum ++;

        mCycleTotalBase[i]++;
        mCycleTotalQual[i]+=qual;

    }

    mReads++;
}

int Stats::getCycles() {
    return mCycles;
}

void Stats::print() {
    if(!summarized) {
        summarize();
    }
    cout << "total reads: " << mReads << endl;
    cout << "total bases: " << mBases << endl;
    cout << "Q20 bases: " << mQ20Total << "(" << (mQ20Total*100.0)/mBases << "%)" << endl;
    cout << "Q30 bases: " << mQ30Total << "(" << (mQ30Total*100.0)/mBases << "%)" << endl;
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
                s->mCycleQ30Bases[i][j] += list[t]->mCycleQ30Bases[i][j];
                s->mCycleQ20Bases[i][j] += list[t]->mCycleQ20Bases[i][j];
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

