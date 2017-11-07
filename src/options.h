#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

class SplitOptions {
public:
    SplitOptions() {
        enabled = false;
        number = 0;
        size = 0;
        digits = 4;
    }
public:
    bool enabled;
    // number of files
    int number;
    // how many reads per file
    int size;
    // digits number of file name prefix, for example 0001 means 4 digits
    int digits;
};

class AdapterOptions {
public:
    AdapterOptions() {
        enabled = true;
    }
public:
    bool enabled;
};

class TrimmingOptions {
public:
    TrimmingOptions() {
        front1 = 0;
        tail1 = 0;
        front2 = 0;
        tail2 = 0;
    }
public:
    // trimming first cycles for read1
    int front1;
    // trimming last cycles for read1
    int tail1;
    // trimming first cycles for read2
    int front2;
    // trimming last cycles for read2
    int tail2;
};

class QualityFilteringOptions {
public:
    QualityFilteringOptions() {
        enabled = true;
        // '0' = Q15
        qualifiedQual = '0';
        unqualifiedPercentLimit = 40;
        nBaseLimit = 5;
    }
public:
    // quality filter enabled
    bool enabled;
    // if a base's quality phred score < qualifiedPhred, then it's considered as a low_qual_base
    char qualifiedQual;
    // if low_qual_base_num > lowQualLimit, then discard this read
    int unqualifiedPercentLimit;
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
    bool validate();
    bool adapterCuttingEnabled();

public:
    // file name of read1 input
    string in1;
    // file name of read2 input
    string in2;
    // file name of read1 output
    string out1;
    // file name of read1 output
    string out2;
    // json file
    string jsonFile;
    // html file
    string htmlFile;
    // compression level
    int compression;
    // worker thread number
    int thread;
    // trimming options
    TrimmingOptions trim;
    // quality filtering options
    QualityFilteringOptions qualfilter;
    // length filtering options
    ReadLengthFilteringOptions lengthFilter;
    // adapter options
    AdapterOptions adapter;
    // multiple file splitting options
    SplitOptions split;

};

#endif