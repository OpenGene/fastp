#ifndef ADAPTER_TRIMMER_H
#define ADAPTER_TRIMMER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "overlapanalysis.h"
#include "filterresult.h"
#include "options.h"

using namespace std;

class AdapterTrimmer{
public:
    AdapterTrimmer();
    ~AdapterTrimmer();

    static bool trimByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr, int diffLimit, int overlapRequire, double diffPercentLimit);
    static bool trimByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr, OverlapResult ov, int frontTrimmed1 = 0, int frontTrimmed2 = 0);
    static bool trimBySequence(Read* r1, FilterResult* fr, string& adapter, bool isR2 = false, int matchReq = 4);
    static bool trimByMultiSequences(Read* r1, FilterResult* fr, vector<string>& adapterList, bool isR2 = false, bool incTrimmedCounter = true);
    static bool test();


};


#endif