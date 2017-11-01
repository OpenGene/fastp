#ifndef STATS_H
#define STATS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

class Stats{
public:
    // this @guessedCycles parameter should be calculated using the first several records
    Stats(int guessedCycles, int bufferMargin = 1000);
    ~Stats();
    int getCycles();

    static Stats* merge(vector<Stats*>& list);

public:
    long mReads;
    /* 
    why we use 8 here?
    map A/T/C/G/N to 0~7 by their ASCII % 8:
    'A' % 8 = 1
    'T' % 8 = 4
    'C' % 8 = 3
    'G' % 8 = 7
    'N' % 8 = 6
    */
    long *mQ30Bases[8];
    long *mQ20Bases[8];
    long *mCycleBaseContents[8];
    long *mCycleTotalBase;
    long *mCycleTotalQual;

private:
    int mCycles;
};

#endif