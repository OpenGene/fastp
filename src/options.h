#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>

using namespace std;

#define UMI_LOC_NONE 0
#define UMI_LOC_INDEX1 1
#define UMI_LOC_INDEX2 2
#define UMI_LOC_READ1 3
#define UMI_LOC_READ2 4
#define UMI_LOC_PER_INDEX 5
#define UMI_LOC_PER_READ 6

class DuplicationOptions {
public:
    DuplicationOptions() {
        enabled = true;
        keylen = 12;
        histSize = 32;
    }
public:
    bool enabled;
    int keylen;
    int histSize;
};

class IndexFilterOptions {
public:
    IndexFilterOptions() {
        enabled = false;
        threshold = 0;
    }
public:
    vector<string> blacklist1;
    vector<string> blacklist2;
    bool enabled;
    int threshold;
};

class LowComplexityFilterOptions {
public:
    LowComplexityFilterOptions() {
        enabled = false;
        threshold = 0.3;
    }
public:
    bool enabled;
    double threshold;
};

class OverrepresentedSequenceAnasysOptions {
public:
    OverrepresentedSequenceAnasysOptions() {
        enabled = false;
        sampling = 20;
    }
public:
    bool enabled;
    int sampling;
};

class PolyGTrimmerOptions {
public:
    PolyGTrimmerOptions() {
        enabled = false;
        minLen = 10;
    }
public:
    bool enabled;
    int minLen;
};

class PolyXTrimmerOptions {
public:
    PolyXTrimmerOptions() {
        enabled = false;
        minLen = 10;
    }
public:
    bool enabled;
    int minLen;
};

class UMIOptions {
public:
    UMIOptions() {
        enabled = false;
        location = UMI_LOC_NONE;
        length = 0;
        skip = 0;
    }
public:
    bool enabled;
    int location;
    int length;
    int skip;
    string prefix;
    string separator;
};

class CorrectionOptions {
public:
    CorrectionOptions() {
        enabled = false;
    }
public:
    bool enabled;
};

class QualityCutOptions {
public:
    QualityCutOptions() {
        enabled5 = false;
        enabled3 = false;
        windowSize = 4;
        quality = 20;
    }
public:
    // enable 5' cutting by quality
    bool enabled5;
    // enable 3' cutting by quality
    bool enabled3;
    // the sliding window size for cutting by quality
    int windowSize;
    // the mean quality requirement for cutting by quality
    int quality;
};

class SplitOptions {
public:
    SplitOptions() {
        enabled = false;
        needEvaluation = false;
        number = 0;
        size = 0;
        digits = 4;
        byFileNumber = false;
        byFileLines = false;
    }
public:
    bool enabled;
    // number of files
    int number;
    // lines of each file
    long size;
    // digits number of file name prefix, for example 0001 means 4 digits
    int digits;
    // need evaluation?
    bool needEvaluation;
    bool byFileNumber;
    bool byFileLines;
};

class AdapterOptions {
public:
    AdapterOptions() {
        enabled = true;
        hasSeqR1 = false;
        hasSeqR2 = false;
    }
public:
    bool enabled;
    string sequence;
    string sequenceR2;
    string detectedAdapter1;
    string detectedAdapter2;
    bool hasSeqR1;
    bool hasSeqR2;
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
        requiredLength = 15;
        maxLength = 0;
    }
public:
    // length filter enabled
    bool enabled;
    // if read_length < requiredLength, then this read is discard
    int requiredLength;
    // length limit, 0 for no limitation
    int maxLength;
};

class Options{
public:
    Options();
    void init();
    bool isPaired();
    bool validate();
    bool adapterCuttingEnabled();
    string getAdapter1();
    string getAdapter2();
    void initIndexFiltering(string blacklistFile1, string blacklistFile2, int threshold = 0);
    vector<string> makeListFromFileByLine(string filename);

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
    // html report title
    string reportTitle;
    // compression level
    int compression;
    // the input file is using phred64 quality scoring
    bool phred64;
    // do not rewrite existing files
    bool dontOverwrite;
    // write STDOUT
    bool outputToSTDOUT;
    // the input R1 file is interleaved
    bool interleavedInput;
    // only process first N reads
    int readsToProcess;
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
    // options for quality cutting
    QualityCutOptions qualityCut;
    // options for base correction
    CorrectionOptions correction;
    // options for UMI
    UMIOptions umi;
    // 3' end polyG trimming, default for Illumina NextSeq/NovaSeq
    PolyGTrimmerOptions polyGTrim;
    // 3' end polyX trimming
    PolyXTrimmerOptions polyXTrim;
    // for overrepresentation analysis
    OverrepresentedSequenceAnasysOptions overRepAnalysis;
    map<string, long> overRepSeqs1;
    map<string, long> overRepSeqs2;
    int seqLen1;
    int seqLen2;
    // low complexity filtering
    LowComplexityFilterOptions complexityFilter;
    // black lists for filtering by index
    IndexFilterOptions indexFilter;
    // options for duplication profiling
    DuplicationOptions duplicate;
    // max value of insert size
    int insertSizeMax;
    // overlap analysis threshold
    int overlapRequire;
    int overlapDiffLimit;

};

#endif