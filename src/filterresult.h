#ifndef FILTER_RESULT_H
#define FILTER_RESULT_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "common.h"
#include "options.h"
#include <fstream>

using namespace std;

class FilterResult{
public:
    FilterResult(Options* opt, bool paired = false);
    ~FilterResult();
    inline long* getFilterReadStats() {return mFilterReadStats;}
    void addFilterResult(int result);
    static FilterResult* merge(vector<FilterResult*>& list);
    void print();
    // for single end
    void addAdapterTrimmed(string adapter);
    // for paired end
    void addAdapterTrimmed(string adapter1, string adapter2);
    // a port of JSON report
    void reportJson(ofstream& ofs, string padding);

public:
    Options* mOptions;
    bool mPaired;
private:
    long mFilterReadStats[FILTER_RESULT_TYPES];
    long mTrimmedAdapterRead;
    long mTrimmedAdapterBases;
};

#endif