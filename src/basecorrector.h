#ifndef BASE_CORRECTOR_H
#define BASE_CORRECTOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "overlapanalysis.h"
#include "filterresult.h"
#include "options.h"

using namespace std;

class BaseCorrector{
public:
    BaseCorrector();
    ~BaseCorrector();

    static int correctByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr, int diffLimit, int overlapRequire, double diffPercentLimit);
    static int correctByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr, OverlapResult ov);
    static bool test();
};


#endif