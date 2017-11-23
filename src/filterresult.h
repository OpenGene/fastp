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
#include <map>

struct classcomp {
    bool operator() (const string& lhs, const string& rhs) const {
        if (lhs.length() < rhs.length())
            return true;
        else if(lhs.length() == rhs.length()) {
            return lhs < rhs;
        } else
            return false;
    }
};

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
    // a part of JSON report
    void reportJson(ofstream& ofs, string padding);
    // a part of JSON report for adapters
    void reportAdapterJson(ofstream& ofs, string padding);
    // a part of HTML report
    void reportHtml(ofstream& ofs, long totalReads, long totalBases);
    // a part of HTML report for adapters
    void reportAdapterHtml(ofstream& ofs, long totalBases);
    void outputAdaptersJson(ofstream& ofs, map<string, long, classcomp>& adapterCounts);
    void outputAdaptersHtml(ofstream& ofs, map<string, long, classcomp>& adapterCounts, long totalBases);
    // deal with base correction results
    long* getCorrectionMatrix() {return mCorrectionMatrix;}
    long getTotalCorrectedBases();
    void addCorrection(char from, char to);
    long getCorrectionNum(char from, char to);
    void incCorrectedReads(int count);

public:
    Options* mOptions;
    bool mPaired;
    long mCorrectedReads;
private:
    long mFilterReadStats[FILTER_RESULT_TYPES];
    long mTrimmedAdapterRead;
    long mTrimmedAdapterBases;
    map<string, long, classcomp> mAdapter1;
    map<string, long, classcomp> mAdapter2;
    long* mCorrectionMatrix;
};

#endif