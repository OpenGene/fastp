#include "options.h"
#include "util.h"
#include <iostream>
#include <fstream>
#include <string.h>

Options::Options(){
    in1 = "";
    in2 = "";
    out1 = "";
    out2 = "";
    reportTitle = "fastp report";
    thread = 1;
    compression = 2;
    phred64 = false;
    dontOverwrite = false;
    inputFromSTDIN = false;
    outputToSTDOUT = false;
    readsToProcess = 0;
    interleavedInput = false;
    insertSizeMax = 512;
    overlapRequire = 30;
    overlapDiffLimit = 5;
    verbose = false;
    seqLen1 = 151;
    seqLen2 = 151;
}

void Options::init() {
}

bool Options::isPaired() {
    return in2.length() > 0 || interleavedInput;
}

bool Options::adapterCuttingEnabled() {
    if(adapter.enabled){
        if(isPaired() || !adapter.sequence.empty())
            return true;
    }
    return false;
}

bool Options::polyXTrimmingEnabled() {
    return polyXTrim.enabled;
}

bool Options::validate() {
    if(in1.empty()) {
        if(!in2.empty())
            error_exit("read2 input is specified by <in2>, but read1 input is not specified by <in1>");
        if(inputFromSTDIN)
            in1 = "/dev/stdin";
        else
            error_exit("read1 input should be specified by --in1, or enable --stdin if you want to read STDIN");
    } else {
        check_file_valid(in1);
    }

    if(!in2.empty()) {
        check_file_valid(in2);
    }

    if(merge.enabled) {
        if(split.enabled) {
            error_exit("splitting mode cannot work with merging mode");
        }
        if(in2.empty() && !interleavedInput) {
            error_exit("read2 input should be specified by --in2 for merging mode");
        }
        if(!out2.empty()) {
            cerr << "Merging mode. Ignore argument --out2 = " << out2 << endl;
            out2 = "";
        }
    }

    // if output to STDOUT, then...
    if(outputToSTDOUT) {
        cerr << "Streaming uncompressed output to STDOUT..." << endl;
        if(!in1.empty() && !in2.empty())
            cerr << "Enable interleaved output mode for paired-end input." << endl;
        if(!out1.empty()) {
            cerr << "Ignore argument --out1 = " << out1 << endl;
            out1 = "";
        }
        if(!out2.empty()) {
            cerr << "Ignore argument --out2 = " << out2 << endl;
            out2 = "";
        }
        if(split.enabled) {
            cerr << "Ignore split mode" << endl;
            split.enabled = false;
        }
        cerr << endl;
    }

    if(in2.empty() && !interleavedInput && !out2.empty()) {
        error_exit("read2 output is specified (--out2), but neighter read2 input is not specified (--in2), nor read1 is interleaved.");
    }

    if(!in2.empty() || interleavedInput) {
        if(!out1.empty() && out2.empty() && !merge.enabled) {
            error_exit("paired-end input, read1 output should be specified together with read2 output (--out2 needed) ");
        }
        if(out1.empty() && !out2.empty()) {
            if(!merge.enabled)
                error_exit("paired-end input, read1 output should be specified (--out1 needed) together with read2 output ");
        }
    }

    if(!in2.empty() && interleavedInput) {
        error_exit("<in2> is not allowed when <in1> is specified as interleaved mode by (--interleaved_in)");
    }

    if(!out1.empty()) {
        //check_file_writable(out1);
        if(out1 == out2) {
            error_exit("read1 output (--out1) and read1 output (--out2) should be different");
        }
        if(dontOverwrite && file_exists(out1)) {
            error_exit(out1 + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
    }
    if(!out2.empty()) {
        //check_file_writable(out2);
        if(dontOverwrite && file_exists(out2)) {
            error_exit(out2 + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
    }

    if(dontOverwrite) {
        if(file_exists(jsonFile)) {
            error_exit(jsonFile + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
        if(file_exists(htmlFile)) {
            error_exit(htmlFile + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
    }

    if(compression < 1 || compression > 9)
        error_exit("compression level (--compression) should be between 1 ~ 9, 1 for fastest, 9 for smallest");

    if(readsToProcess < 0)
        error_exit("the number of reads to process (--reads_to_process) cannot be negative");

    if(thread < 1) {
        thread = 1;
    } else if(thread > 16) {
        cerr << "WARNING: fastp uses up to 16 threads although you specified " << thread << endl;
        thread = 16;
    }

    if(trim.front1 < 0 || trim.front1 > 30)
        error_exit("trim_front1 (--trim_front1) should be 0 ~ 30, suggest 0 ~ 4");

    if(trim.tail1 < 0 || trim.tail1 > 100)
        error_exit("trim_tail1 (--trim_tail1) should be 0 ~ 100, suggest 0 ~ 4");

    if(trim.front2 < 0 || trim.front2 > 30)
        error_exit("trim_front2 (--trim_front2) should be 0 ~ 30, suggest 0 ~ 4");

    if(trim.tail2 < 0 || trim.tail2 > 100)
        error_exit("trim_tail2 (--trim_tail2) should be 0 ~ 100, suggest 0 ~ 4");

    if(qualfilter.qualifiedQual - 33 < 0 || qualfilter.qualifiedQual - 33 > 93)
        error_exit("qualitified phred (--qualified_quality_phred) should be 0 ~ 93, suggest 10 ~ 20");

    if(qualfilter.unqualifiedPercentLimit < 0 || qualfilter.unqualifiedPercentLimit > 100)
        error_exit("unqualified percent limit (--unqualified_percent_limit) should be 0 ~ 100, suggest 20 ~ 60");

    if(qualfilter.nBaseLimit < 0 || qualfilter.nBaseLimit > 50)
        error_exit("N base limit (--n_base_limit) should be 0 ~ 50, suggest 3 ~ 10");

    if(lengthFilter.requiredLength < 0 )
        error_exit("length requirement (--length_required) should be >0, suggest 15 ~ 100");

    if(split.enabled ) {
        if(split.digits < 0 || split.digits > 10)
            error_exit("you have enabled splitting output to multiple files, the digits number of file name prefix (--split_prefix_digits) should be 0 ~ 10.");
        
        if(split.byFileNumber) {
            if(split.number < 2 || split.number >= 1000)
                error_exit("you have enabled splitting output by file number, the number of files (--split) should be 2 ~ 999.");
            // thread number cannot be more than the number of file to split
            if(thread > split.number)
                thread = split.number;
        }

        if(split.byFileLines) {
            if(split.size < 1000/4)
                error_exit("you have enabled splitting output by file lines, the file lines (--split_by_lines) should be >= 1000.");
        }
    }

    if(qualityCut.enabledFront || qualityCut.enabledTail || qualityCut.enabledRight) {
        if(qualityCut.windowSizeShared < 1 || qualityCut.windowSizeShared > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_window_size) should be between 1~1000.");
        if(qualityCut.qualityShared < 1 || qualityCut.qualityShared > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_mean_quality) should be 1 ~ 30, suggest 15 ~ 20.");
        if(qualityCut.windowSizeFront < 1 || qualityCut.windowSizeFront > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_front_window_size) should be between 1~1000.");
        if(qualityCut.qualityFront < 1 || qualityCut.qualityFront > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_front_mean_quality) should be 1 ~ 30, suggest 15 ~ 20.");
        if(qualityCut.windowSizeTail < 1 || qualityCut.windowSizeTail > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_tail_window_size) should be between 1~1000.");
        if(qualityCut.qualityTail < 1 || qualityCut.qualityTail > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_tail_mean_quality) should be 1 ~ 30, suggest 13 ~ 20.");
        if(qualityCut.windowSizeRight < 1 || qualityCut.windowSizeRight > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_right_window_size) should be between 1~1000.");
        if(qualityCut.qualityRight < 1 || qualityCut.qualityRight > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_right_mean_quality) should be 1 ~ 30, suggest 15 ~ 20.");
    }

    if(adapter.sequence!="auto" && !adapter.sequence.empty()) {
        // validate adapter sequence for single end adapter trimming
        if(adapter.sequence.length() <= 3)
            error_exit("the sequence of <adapter_sequence> should be longer than 3");

        // validate bases
        for(int i=0; i<adapter.sequence.length(); i++) {
            char c = adapter.sequence[i];
            if(c!='A' && c!='T' && c!='C' && c!='G') {
                error_exit("the adapter <adapter_sequence> can only have bases in {A, T, C, G}, but the given sequence is: " + adapter.sequence);
            }
        }

        adapter.hasSeqR1 = true;
    }

    if(adapter.sequenceR2!="auto" && !adapter.sequenceR2.empty()) {
        // validate adapter sequenceR2 for single end adapter trimming
        if(adapter.sequenceR2.length() <= 3)
            error_exit("the sequence of <adapter_sequence_r2> should be longer than 3");

        // validate bases
        for(int i=0; i<adapter.sequenceR2.length(); i++) {
            char c = adapter.sequenceR2[i];
            if(c!='A' && c!='T' && c!='C' && c!='G') {
                error_exit("the adapter <adapter_sequence_r2> can only have bases in {A, T, C, G}, but the given sequenceR2 is: " + adapter.sequenceR2);
            }
        }

        adapter.hasSeqR2 = true;
    }

    if(correction.enabled && !isPaired()) {
        cerr << "WARNING: base correction is only appliable for paired end data, ignored -c/--correction" << endl;
        correction.enabled = false;
    }

    if(umi.enabled) {
        if(umi.location == UMI_LOC_READ1 || umi.location == UMI_LOC_READ2 || umi.location == UMI_LOC_PER_READ) {
            if(umi.length<1 || umi.length>100)
                error_exit("UMI length should be 1~100");
            if(umi.skip<0 || umi.skip>100)
                error_exit("The base number to skip after UMI <umi_skip> should be 0~100");
        }else {
            if(umi.skip>0)
                error_exit("Only if the UMI location is in read1/read2/per_read, you can skip bases after UMI");
            if(umi.length>0)
                error_exit("Only if the UMI location is in read1/read2/per_read, you can set the UMI length");
        }
        if(!umi.prefix.empty()) {
            if(umi.prefix.length() >= 10)
                error_exit("UMI prefix should be shorter than 10");
            for(int i=0; i<umi.prefix.length(); i++) {
                char c = umi.prefix[i];
                if( !(c>='A' && c<='Z') && !(c>='a' && c<='z') && !(c>='0' && c<='9')) {
                    error_exit("UMI prefix can only have characters and numbers, but the given is: " + umi.prefix);
                }
            }
        }
        if(!umi.separator.empty()) {
            if(umi.separator.length()>10)
                error_exit("UMI separator cannot be longer than 10 base pairs");
            // validate bases
            for(int i=0; i<umi.separator.length(); i++) {
                char c = umi.separator[i];
                if(c!='A' && c!='T' && c!='C' && c!='G') {
                    error_exit("UMI separator can only have bases in {A, T, C, G}, but the given sequence is: " + umi.separator);
                }
            }
        }

    }

    if(overRepAnalysis.sampling < 1 || overRepAnalysis.sampling > 10000)
        error_exit("overrepresentation_sampling should be 1~10000");

    return true;
}

bool Options::shallDetectAdapter(bool isR2) {
    if(!adapter.enabled)
        return false;

    if(isR2) {
        return isPaired() && adapter.detectAdapterForPE && adapter.sequenceR2 == "auto";
    } else {
        if(isPaired())
            return adapter.detectAdapterForPE && adapter.sequence == "auto";
        else
            return adapter.sequence == "auto";
    }
}

void Options::initIndexFiltering(string blacklistFile1, string blacklistFile2, int threshold) {
    if(blacklistFile1.empty() && blacklistFile2.empty())
        return;

    if(!blacklistFile1.empty()){
        check_file_valid(blacklistFile1);
        indexFilter.blacklist1 = makeListFromFileByLine(blacklistFile1);
    }

    if(!blacklistFile2.empty()){
        check_file_valid(blacklistFile2);
        indexFilter.blacklist2 = makeListFromFileByLine(blacklistFile2);
    }

    if(indexFilter.blacklist1.empty() && indexFilter.blacklist2.empty())
        return;

    indexFilter.enabled = true;
    indexFilter.threshold = threshold;
}

vector<string> Options::makeListFromFileByLine(string filename) {
    vector<string> ret;
    ifstream file;
    file.open(filename.c_str(), ifstream::in);
    const int maxLine = 1000;
    char line[maxLine];
    cerr << "filter by index, loading " << filename << endl;
    while(file.getline(line, maxLine)){
        // trim \n, \r or \r\n in the tail
        int readed = strlen(line);
        if(readed >=2 ){
            if(line[readed-1] == '\n' || line[readed-1] == '\r'){
                line[readed-1] = '\0';
                if(line[readed-2] == '\r')
                    line[readed-2] = '\0';
            }
        }
        string linestr(line);
        for(int i=0; i<linestr.length(); i++) {
            if(linestr[i] != 'A' && linestr[i] != 'T' && linestr[i] != 'C' && linestr[i] != 'G') {
                error_exit("processing " + filename + ", each line should be one barcode, which can only contain A/T/C/G");
            }
        }
        cerr << linestr << endl;
        ret.push_back(linestr);
    }
    cerr << endl;
    return ret;
}

string Options::getAdapter1(){
    if(adapter.sequence == "" || adapter.sequence == "auto")
        return "unspecified";
    else
        return adapter.sequence;
}

string Options::getAdapter2(){
    if(adapter.sequenceR2 == "" || adapter.sequenceR2 == "auto")
        return "unspecified";
    else
        return adapter.sequenceR2;
}
