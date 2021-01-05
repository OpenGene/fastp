#ifndef OVERLAP_ANALYSIS_H
#define OVERLAP_ANALYSIS_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "common.h"
#include "options.h"
#include "read.h"

using namespace std;

class OverlapResult {
public:
    bool overlapped;
    int offset;
    int overlap_len;
    int diff;
};

class OverlapAnalysis{
public:
    OverlapAnalysis();
    ~OverlapAnalysis();

    static OverlapResult analyze(Sequence&  r1, Sequence&  r2, int diffLimit, int overlapRequire, double diffPercentLimit);
    static OverlapResult analyze(Read* r1, Read* r2, int diffLimit, int overlapRequire, double diffPercentLimit);
    static Read* merge(Read* r1, Read* r2, OverlapResult ov);

public:
    static bool test();

};

#endif