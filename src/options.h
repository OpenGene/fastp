#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

class TrimmingOptions {
public:
    TrimmingOptions() {
        head = 0;
        tail = 0;
    }
public:
    // trimming first cycles
    int head;
    // trimming last cycles
    int tail;
};

class QualityFilteringOptions {
public:
    QualityFilteringOptions() {
        qualifiedPhred = 15;
        unqualifiedBaseLimit = 5;
        nBaseLimit = 3;
    }
public:
    // if a base's quality phred score < qualifiedPhred, then it's considered as a low_qual_base
    int qualifiedPhred;
    // if low_qual_base_num > lowQualLimit, then discard this read
    int unqualifiedBaseLimit;
    // if n_base_number > nBaseLimit, then discard this read
    int nBaseLimit;
};

class ReadLengthFilteringOptions {
public:
    ReadLengthFilteringOptions() {
        enabled = false;
        requiredLength = 30;
    }
public:
    // length filter enabled
    bool enabled;
    // if read_length < requiredLength, then this read is discard
    int requiredLength;
};

class Options{
public:
    Options();
    void init();
    bool isPaired();

public:
    // file name of read1 input
    string in1;
    // file name of read2 input
    string in2;
    // file name of read1 output
    string out1;
    // file name of read1 output
    string out2;
    // worker thread number
    int thread;
    // trimming options
    TrimmingOptions trim;
    // quality filtering options
    QualityFilteringOptions qualfilter;

};

#endif