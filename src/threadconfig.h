#ifndef THREAD_CONFIG_H
#define THREAD_CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "stats.h"

using namespace std;

class ThreadConfig{
public:
    ThreadConfig(int seqCycles, bool paired = false);
    inline Stats* getLeftReadStats() {return mStats1;}

public:
    Stats* mStats1;
    Stats* mStats2;
};

#endif