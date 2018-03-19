#include <stdio.h>
#include "fastqreader.h"
#include "unittest.h"
#include <time.h>
#include "cmdline.h"
#include <sstream>
#include "util.h"
#include "options.h"
#include "processor.h"
#include "evaluator.h"

string command;

int main(int argc, char* argv[]){
    // display version info if no argument is given
    if(argc == 1) {
        cout << "fastp: an ultra-fast all-in-one FASTQ preprocessor" << endl << "version " << FASTP_VER << endl;
    }
    if (argc == 2 && strcmp(argv[1], "test")==0){
        UnitTest tester;
        return tester.run();
    }
    if (argc == 2 && (strcmp(argv[1], "-v")==0 || strcmp(argv[1], "--version")==0)){
        cout << "fastp: an ultra-fast all-in-one FASTQ preprocessor" << endl << "version " << FASTP_VER << endl;
        return 0;
    }
    cmdline::parser cmd;
    // input/output
    cmd.add<string>("in1", 'i', "read1 input file name", true, "");
    cmd.add<string>("out1", 'o', "read1 output file name", false, "");
    cmd.add<string>("in2", 'I', "read2 input file name", false, "");
    cmd.add<string>("out2", 'O', "read2 output file name", false, "");
    cmd.add("phred64", '6', "indicates the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)");
    cmd.add<int>("compression", 'z', "compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 2.", false, 2);

    // adapter
    cmd.add("disable_adapter_trimming", 'A', "adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled");
    cmd.add<string>("adapter_sequence", 'a', "the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped.", false, "auto");
    cmd.add<string>("adapter_sequence_r2", 0, "the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence>", false, "");

    // trimming
    cmd.add<int>("trim_front1", 'f', "trimming how many bases in front for read1, default is 0", false, 0);
    cmd.add<int>("trim_tail1", 't', "trimming how many bases in tail for read1, default is 0", false, 0);
    cmd.add<int>("trim_front2", 'F', "trimming how many bases in front for read2. If it's not specified, it will follow read1's settings", false, 0);
    cmd.add<int>("trim_tail2", 'T', "trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings", false, 0);

    // polyG tail trimming
    cmd.add("trim_poly_g", 'g', "force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data");
    cmd.add<int>("poly_g_min_len", 0, "the minimum length to detect polyG in the read tail. 10 by default.", false, 10);
    cmd.add("disable_trim_poly_g", 'G', "disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data");
    
    // polyX tail trimming
    cmd.add("trim_poly_x", 'x', "enable polyX trimming in 3' ends.");
    cmd.add<int>("poly_x_min_len", 0, "the minimum length to detect polyX in the read tail. 10 by default.", false, 10);

    // sliding window cutting for each reads
    cmd.add("cut_by_quality5", '5', "enable per read cutting by quality in front (5'), default is disabled (WARNING: this will interfere deduplication for both PE/SE data)");
    cmd.add("cut_by_quality3", '3', "enable per read cutting by quality in tail (3'), default is disabled (WARNING: this will interfere deduplication for SE data)");
    cmd.add<int>("cut_window_size", 'W', "the size of the sliding window for sliding window trimming, default is 4", false, 4);
    cmd.add<int>("cut_mean_quality", 'M', "the bases in the sliding window with mean quality below cutting_quality will be cut, default is Q20", false, 20);

    // quality filtering
    cmd.add("disable_quality_filtering", 'Q', "quality filtering is enabled by default. If this option is specified, quality filtering is disabled");
    cmd.add<int>("qualified_quality_phred", 'q', "the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified.", false, 15);
    cmd.add<int>("unqualified_percent_limit", 'u', "how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%", false, 40);
    cmd.add<int>("n_base_limit", 'n', "if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5", false, 5);

    // length filtering
    cmd.add("disable_length_filtering", 'L', "length filtering is enabled by default. If this option is specified, length filtering is disabled");
    cmd.add<int>("length_required", 'l', "reads shorter than length_required will be discarded, default is 15.", false, 15);

    // low complexity filtering
    cmd.add("low_complexity_filter", 'y', "enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).");
    cmd.add<int>("complexity_threshold", 'Y', "the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required.", false, 30);
    
    // base correction in overlapped regions of paired end data
    cmd.add("correction", 'c', "enable base correction in overlapped regions (only for PE data), default is disabled");

    // umi
    cmd.add("umi", 'U', "enable unique molecular identifer (UMI) preprocessing");
    cmd.add<string>("umi_loc", 0, "specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none", false, "");
    cmd.add<int>("umi_len", 0, "if the UMI is in read1/read2, its length should be provided", false, 0);
    cmd.add<string>("umi_prefix", 0, "if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default", false, "");
    //cmd.add<string>("umi_sep", 0, "if the UMI is in read1/read2, it can have a separator (several bases separate UMI and insert DNA), default is empty", false, "");

    // overrepresented sequence analysis
    cmd.add("overrepresentation_analysis", 'p', "enable overrepresented sequence analysis.");
    cmd.add<int>("overrepresentation_sampling", 'P', "one in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis (1~10000), smaller is slower, default is 20.", false, 20);
    
    // reporting
    cmd.add<string>("json", 'j', "the json format report file name", false, "fastp.json");
    cmd.add<string>("html", 'h', "the html format report file name", false, "fastp.html");
    cmd.add<string>("report_title", 'R', "should be quoted with \' or \", default is \"fastp report\"", false, "fastp report");

    // threading
    cmd.add<int>("thread", 'w', "worker thread number, default is 3", false, 3);

    // split the output
    cmd.add<int>("split", 's', "split output by limiting total split file number with this option (2~999), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default", false, 0);
    cmd.add<long>("split_by_lines", 'S', "split output by limiting lines of each file with this option(>=1000), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default", false, 0);
    cmd.add<int>("split_prefix_digits", 'd', "the digits for the sequential number padding (1~10), default is 4, so the filename will be padded as 0001.xxx, 0 to disable padding", false, 4);

    cmd.parse_check(argc, argv);

    Options opt;

    // I/O
    opt.in1 = cmd.get<string>("in1");
    opt.in2 = cmd.get<string>("in2");
    opt.out1 = cmd.get<string>("out1");
    opt.out2 = cmd.get<string>("out2");
    opt.compression = cmd.get<int>("compression");
    opt.phred64 = cmd.exist("phred64");

    // adapter cutting
    opt.adapter.enabled = !cmd.exist("disable_adapter_trimming");
    opt.adapter.sequence = cmd.get<string>("adapter_sequence");
    opt.adapter.sequenceR2 = cmd.get<string>("adapter_sequence_r2");
    if(opt.adapter.sequenceR2.empty() && opt.adapter.sequence != "auto") {
        opt.adapter.sequenceR2 = opt.adapter.sequence;
    }

    // trimming
    opt.trim.front1 = cmd.get<int>("trim_front1");
    opt.trim.tail1 = cmd.get<int>("trim_tail1");
    // read2 settings follows read1 if it's not specified
    if(cmd.exist("trim_front2"))
        opt.trim.front2 = cmd.get<int>("trim_front2");
    else
        opt.trim.front2 = opt.trim.front1;
    if(cmd.exist("trim_tail2"))
        opt.trim.tail2 = cmd.get<int>("trim_tail2");
    else
        opt.trim.tail2 = opt.trim.tail1;

    // polyG tail trimming
    if(cmd.exist("trim_poly_g") && cmd.exist("disable_trim_poly_g")) {
        error_exit("You cannot enabled both trim_poly_g and disable_trim_poly_g");
    } else if(cmd.exist("trim_poly_g")) {
        opt.polyGTrim.enabled = true;
    } else if(cmd.exist("disable_trim_poly_g")) {
        opt.polyGTrim.enabled = false;
    }
    opt.polyGTrim.minLen = cmd.get<int>("poly_g_min_len");

    // polyX tail trimming
    if(cmd.exist("trim_poly_x")) {
        opt.polyXTrim.enabled = true;
    }
    opt.polyXTrim.minLen = cmd.get<int>("poly_x_min_len");

    // sliding window cutting by quality
    opt.qualityCut.enabled5 = cmd.exist("cut_by_quality5");
    opt.qualityCut.enabled3 = cmd.exist("cut_by_quality3");
    opt.qualityCut.windowSize = cmd.get<int>("cut_window_size");
    opt.qualityCut.quality = cmd.get<int>("cut_mean_quality");
    // raise a warning if -5/-3 is not enabled but -W/-M is enabled
    if(cmd.exist("cut_window_size") && !opt.qualityCut.enabled5 && !opt.qualityCut.enabled3) {
        cerr << "WARNING: you've specified sliding window size (-W/--cut_window_size), but you haven't enabled per read cutting by quality for 5'(-5/--cut_by_quality5) or 3' (-3/--cut_by_quality3), so quality cutting is ignored " << endl << endl;
    } else if(cmd.exist("cut_mean_quality") && !opt.qualityCut.enabled5 && !opt.qualityCut.enabled3) {
        cerr << "WARNING: you've specified sliding window mean quality requirement (-M/--cut_mean_quality), but you haven't enabled per read cutting by quality for 5'(-5/--cut_by_quality5) or 3' (-3/--cut_by_quality3), so quality cutting is ignored "<<endl << endl;
    }

    // quality filtering
    opt.qualfilter.enabled = !cmd.exist("disable_quality_filtering");
    opt.qualfilter.qualifiedQual = num2qual(cmd.get<int>("qualified_quality_phred"));
    opt.qualfilter.unqualifiedPercentLimit = cmd.get<int>("unqualified_percent_limit");
    opt.qualfilter.nBaseLimit = cmd.get<int>("n_base_limit");

    // length filtering
    opt.lengthFilter.enabled = !cmd.exist("disable_length_filtering");
    opt.lengthFilter.requiredLength = cmd.get<int>("length_required");

    // low complexity filter
    opt.complexityFilter.enabled = cmd.exist("low_complexity_filter");
    opt.complexityFilter.threshold = (min(100, max(0, cmd.get<int>("complexity_threshold")))) / 100.0;

    // overlap correction
    opt.correction.enabled = cmd.exist("correction");

    // threading
    opt.thread = cmd.get<int>("thread");

    // reporting
    opt.jsonFile = cmd.get<string>("json");
    opt.htmlFile = cmd.get<string>("html");
    opt.reportTitle = cmd.get<string>("report_title");

    // splitting
    opt.split.enabled = cmd.exist("split") || cmd.exist("split_by_lines");
    opt.split.digits = cmd.get<int>("split_prefix_digits");
    if(cmd.exist("split") && cmd.exist("split_by_lines")) {
        error_exit("You cannot set both splitting by file number (--split) and splitting by file lines (--split_by_lines), please choose either.");
    }
    if(cmd.exist("split")) {
        opt.split.number = cmd.get<int>("split");
        opt.split.needEvaluation = true;
        opt.split.byFileNumber = true;
    }
    if(cmd.exist("split_by_lines")) {
        long lines = cmd.get<long>("split_by_lines");
        if(lines % 4 != 0) {
            error_exit("Line number (--split_by_lines) should be a multiple of 4");
        }
        opt.split.size = lines / 4; // 4 lines per record
        opt.split.needEvaluation = false;
        opt.split.byFileLines = true;
    }

    // umi
    opt.umi.enabled = cmd.exist("umi");
    opt.umi.length = cmd.get<int>("umi_len");
    opt.umi.prefix = cmd.get<string>("umi_prefix");
    if(opt.umi.enabled) {
        string umiLoc = cmd.get<string>("umi_loc");
        str2lower(umiLoc);
        if(umiLoc.empty())
            error_exit("You've enabled UMI by (--umi), you should specify the UMI location by (--umi_loc)");
        if(umiLoc != "index1" && umiLoc != "index2" && umiLoc != "read1" && umiLoc != "read2" && umiLoc != "per_index" && umiLoc != "per_read") {
            error_exit("UMI location can only be index1/index2/read1/read2/per_index/per_read");
        }
        if(!opt.isPaired() && (umiLoc == "index2" || umiLoc == "read2"))
            error_exit("You specified the UMI location as " + umiLoc + ", but the input data is not paired end.");
        if(opt.umi.length == 0 && (umiLoc == "read1" || umiLoc == "read2" ||  umiLoc == "per_read"))
            error_exit("You specified the UMI location as " + umiLoc + ", but the length is not specified (--umi_len).");
        if(umiLoc == "index1") {
            opt.umi.location = UMI_LOC_INDEX1;
        } else if(umiLoc == "index2") {
            opt.umi.location = UMI_LOC_INDEX2;
        } else if(umiLoc == "read1") {
            opt.umi.location = UMI_LOC_READ1;
        } else if(umiLoc == "read2") {
            opt.umi.location = UMI_LOC_READ2;
        } else if(umiLoc == "per_index") {
            opt.umi.location = UMI_LOC_PER_INDEX;
        } else if(umiLoc == "per_read") {
            opt.umi.location = UMI_LOC_PER_READ;
        }
    }

    // overrepresented sequence analysis
    opt.overRepAnalysis.enabled = cmd.exist("overrepresentation_analysis");
    opt.overRepAnalysis.sampling = cmd.get<int>("overrepresentation_sampling");

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    time_t t1 = time(NULL);

    Evaluator eva(&opt);
    eva.evaluateSeqLen();

    if(opt.overRepAnalysis.enabled)
        eva.evaluateOverRepSeqs();

    long readNum = 0;

    // using evaluator to guess how many reads in total
    if(opt.adapter.enabled && !opt.isPaired() && opt.adapter.sequence == "auto") {
        cout << "Detecting adapter for single end input..." << endl;
        string adapt = eva.evaluateRead1AdapterAndReadNum(readNum);
        if(adapt.length() >= 12 ) {
            cout << "Detected adapter: " << adapt << endl << endl;
            opt.adapter.sequence = adapt;
        } else {
            cout << "No adapter detected" << endl << endl;
            opt.adapter.sequence = "";
        }
    }

    opt.validate();

    // using evaluator to guess how many reads in total
    if(opt.split.needEvaluation) {
        // if readNum is not 0, means it is already evaluated by other functions
        if(readNum == 0) {
            eva.evaluateReadNum(readNum);
        }
        opt.split.size = readNum / opt.split.number;
        // one record per file at least
        if(opt.split.size <= 0) {
            opt.split.size = 1;
            cerr << "WARNING: the input file has less reads than the number of files to split" << endl;
        }
    }

    // using evaluator to check if it's two color system
    if(!cmd.exist("trim_poly_g") && !cmd.exist("disable_trim_poly_g")) {
        bool twoColorSystem = eva.isTwoColorSystem();
        if(twoColorSystem){
            opt.polyGTrim.enabled = true;
        }
    }

    Processor p(&opt);
    p.process();
    
    time_t t2 = time(NULL);

    cout << endl << "JSON report: " << opt.jsonFile << endl;
    cout << "HTML report: " << opt.htmlFile << endl;
    cout << endl << command << endl;
    cout << "fastp v" << FASTP_VER << ", time used: " << (t2)-t1 << " seconds" << endl;

    return 0;
}
