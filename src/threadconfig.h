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
    inline Stats* getPreStats1() {return mPreStats1;}
    inline Stats* getPostStats1() {return mPostStats1;}
    inline Stats* getPreStats2() {return mPreStats2;}
    inline Stats* getPostStats2() {return mPostStats2;}

public:
    Stats* mPreStats1;
    Stats* mPostStats1;
    Stats* mPreStats2;
    Stats* mPostStats2;
};

#endif