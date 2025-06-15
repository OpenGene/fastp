#ifndef HTML_REPORTER_H
#define HTML_REPORTER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "options.h"
#include "stats.h"
#include "filterresult.h"
#include <fstream>
#include <atomic>
#include "common.h"
#include "util.h"

using namespace std;

class HtmlReporter{
public:
    HtmlReporter(Options* opt);
    ~HtmlReporter();
    void setDup(double dupRate);
    void setInsertHist(atomic_long* insertHist, int insertSizePeak);
    void report(FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2 = NULL, Stats* postStats2 = NULL);

    static void outputRow(ofstream& ofs, string key, long value);
    static void outputRow(ofstream& ofs, string key, string value);
    static string formatNumber(long number);
    static string getPercents(long numerator, long denominator);
private:
    const string getCurrentSystemTime();
    void printHeader(ofstream& ofs);
    void printCSS(ofstream& ofs);
    void printJS(ofstream& ofs);
    void printFooter(ofstream& ofs);
    void reportQualHistogram(ofstream& ofs, string caption, Stats* stats1, Stats* stats2 = NULL);
    void reportDuplication(ofstream& ofs);
    void reportInsertSize(ofstream& ofs, int isizeLimit);
    void printSummary(ofstream& ofs, FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2);
    
private:
    Options* mOptions;
    int* mDupHist;
    double* mDupMeanGC;
    double mDupRate;
    atomic_long* mInsertHist;
    int mInsertSizePeak;
};


#endif
