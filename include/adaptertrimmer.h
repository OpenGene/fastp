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

    static bool trimByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr);
    static bool trimByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr, OverlapResult ov);
    static bool trimBySequence(Read* r1, FilterResult* fr, string& adapter, bool isR2 = false);
    static bool test();


};


#endif