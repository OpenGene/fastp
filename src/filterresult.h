#ifndef FILTER_RESULT_H
#define FILTER_RESULT_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "common.h"
#include "options.h"

using namespace std;

class FilterResult{
public:
    FilterResult(Options* opt, bool paired = false);
    ~FilterResult();
    inline long* getFilterReadStats() {return mFilterReadStats;}
    void addFilterResult(int result);
    static FilterResult* merge(vector<FilterResult*>& list);
    void print();

public:
    Options* mOptions;
    bool mPaired;
private:
    long mFilterReadStats[FILTER_RESULT_TYPES];
};

#endif