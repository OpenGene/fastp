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

// TODO: code refactoring to remove these global variables
string command;
mutex logmtx;

int main(int argc, char* argv[]){
    // display version info if no argument is given
    if(argc == 1) {
        cerr << "fastp: an ultra-fast all-in-one FASTQ preprocessor" << endl << "version " << FASTP_VER << endl;
        //cerr << "fastp --help to see the help"<<endl;
        //return 0;
    }
    if (argc == 2 && strcmp(argv[1], "test")==0){
        UnitTest tester;
        tester.run();
        return 0;
    }
    if (argc == 2 && (strcmp(argv[1], "-v")==0 || strcmp(argv[1], "--version")==0)){
        cout << "fastp " << FASTP_VER << endl;
        return 0;
    }
    cmdline::parser cmd;
    // input/output
    cmd.add<string>("in1", 'i', "read1 input file name", false, "");
    cmd.add<string>("out1", 'o', "read1 output file name", false, "");
    cmd.add<string>("in2", 'I', "read2 input file name", false, "");
    cmd.add<string>("out2", 'O', "read2 output file name", false, "");
    cmd.add<string>("unpaired1", 0, "for PE input, if read1 passed QC but read2 not, it will be written to unpaired1. Default is to discard it.", false, "");
    cmd.add<string>("unpaired2", 0, "for PE input, if read2 passed QC but read1 not, it will be written to unpaired2. If --unpaired2 is same as --unpaired1 (default mode), both unpaired reads will be written to this same file.", false, "");
    cmd.add<string>("overlapped_out", 0, "for each read pair, output the overlapped region if it has no any mismatched base.", false, "");
    cmd.add<string>("failed_out", 0, "specify the file to store reads that cannot pass the filters.", false, "");
    cmd.add("merge", 'm', "for paired-end input, merge each pair of reads into a single read if they are overlapped. The merged reads will be written to the file given by --merged_out, the unmerged reads will be written to the files specified by --out1 and --out2. The merging mode is disabled by default.");
    cmd.add<string>("merged_out", 0, "in the merging mode, specify the file name to store merged output, or specify --stdout to stream the merged output", false, "");
    cmd.add("include_unmerged", 0, "in the merging mode, write the unmerged or unpaired reads to the file specified by --merge. Disabled by default.");
    cmd.add("phred64", '6', "indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)");
    cmd.add<int>("compression", 'z', "compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4.", false, 4);
    cmd.add("stdin", 0, "input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in.");
    cmd.add("stdout", 0, "stream passing-filters reads to STDOUT. This option will result in interleaved FASTQ output for paired-end output. Disabled by default.");
    cmd.add("interleaved_in", 0, "indicate that <in1> is an interleaved FASTQ which contains both read1 and read2. Disabled by default.");
    cmd.add<int>("reads_to_process", 0, "specify how many reads/pairs to be processed. Default 0 means process all reads.", false, 0);
    cmd.add("dont_overwrite", 0, "don't overwrite existing files. Overwritting is allowed by default.");
    cmd.add("fix_mgi_id", 0, "the MGI FASTQ ID format is not compatible with many BAM operation tools, enable this option to fix it.");
    cmd.add("verbose", 'V', "output verbose log information (i.e. when every 1M reads are processed).");

    // adapter
    cmd.add("disable_adapter_trimming", 'A', "adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled");
    cmd.add<string>("adapter_sequence", 'a', "the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped.", false, "auto");
    cmd.add<string>("adapter_sequence_r2", 0, "the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence>", false, "auto");
    cmd.add<string>("adapter_fasta", 0, "specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file", false, "");
    cmd.add("detect_adapter_for_pe", '2', "enable adapter detection for PE data to get ultra-clean data. It takes more time to find just a little bit more adapters.");
    cmd.add("allow_gap_overlap_trimming", 0, "allow up to one gap when trim adapters by overlap analysis for PE data. By default no gap is allowed.");
    cmd.add<int>("dimer_max_len", 0, "if the read length is less than or equal to this value after adapter trimming, it is considered an adapter dimer. Requires adapter evidence.", false, 2);

    // trimming
    cmd.add<int>("trim_front1", 'f', "trimming how many bases in front for read1, default is 0", false, 0);
    cmd.add<int>("trim_tail1", 't', "trimming how many bases in tail for read1, default is 0", false, 0);
    cmd.add<int>("max_len1", 'b', "if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation", false, 0);
    cmd.add<int>("trim_front2", 'F', "trimming how many bases in front for read2. If it's not specified, it will follow read1's settings", false, 0);
    cmd.add<int>("trim_tail2", 'T', "trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings", false, 0);
    cmd.add<int>("max_len2", 'B', "if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings", false, 0);

    // duplication evaluation and deduplication
    cmd.add("dedup", 'D', "enable deduplication to drop the duplicated reads/pairs");
    cmd.add<int>("dup_calc_accuracy", 0, "accuracy level to calculate duplication (1~6), higher level uses more memory (1G, 2G, 4G, 8G, 16G, 32G). Default 1 for no-dedup mode, and 3 for dedup mode.", false);
    cmd.add("dont_eval_duplication", 0, "don't evaluate duplication rate to save time and use less memory.");

    // polyG tail trimming
    cmd.add("trim_poly_g", 'g', "force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data");
    cmd.add<int>("poly_g_min_len", 0, "the minimum length to detect polyG in the read tail. 10 by default.", false, 10);
    cmd.add("disable_trim_poly_g", 'G', "disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data");
    
    // polyX tail trimming
    cmd.add("trim_poly_x", 'x', "enable polyX trimming in 3' ends.");
    cmd.add<int>("poly_x_min_len", 0, "the minimum length to detect polyX in the read tail. 10 by default.", false, 10);

    // cutting by quality
    cmd.add("cut_front", '5', "move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.");
    cmd.add("cut_tail", '3', "move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.");
    cmd.add("cut_right", 'r', "move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.");
    cmd.add<int>("cut_window_size", 'W', "the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4", false, 4);
    cmd.add<int>("cut_mean_quality", 'M', "the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20)", false, 20);
    cmd.add<int>("cut_front_window_size", 0, "the window size option of cut_front, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_front_mean_quality", 0, "the mean quality requirement option for cut_front, default to cut_mean_quality if not specified", false, 20);
    cmd.add<int>("cut_tail_window_size", 0, "the window size option of cut_tail, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_tail_mean_quality", 0, "the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified", false, 20);
    cmd.add<int>("cut_right_window_size", 0, "the window size option of cut_right, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_right_mean_quality", 0, "the mean quality requirement option for cut_right, default to cut_mean_quality if not specified", false, 20);


    // quality filtering
    cmd.add("disable_quality_filtering", 'Q', "quality filtering is enabled by default. If this option is specified, quality filtering is disabled");
    cmd.add<int>("qualified_quality_phred", 'q', "the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified.", false, 15);
    cmd.add<int>("unqualified_percent_limit", 'u', "how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%", false, 40);
    cmd.add<int>("n_base_limit", 'n', "if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5", false, 5);
    cmd.add<int>("average_qual", 'e', "if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement", false, 0);

    // length filtering
    cmd.add("disable_length_filtering", 'L', "length filtering is enabled by default. If this option is specified, length filtering is disabled");
    cmd.add<int>("length_required", 'l', "reads shorter than length_required will be discarded, default is 15.", false, 15);
    cmd.add<int>("length_limit", 0, "reads longer than length_limit will be discarded, default 0 means no limitation.", false, 0);

    // low complexity filtering
    cmd.add("low_complexity_filter", 'y', "enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).");
    cmd.add<int>("complexity_threshold", 'Y', "the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required.", false, 30);

    // filter by indexes
    cmd.add<string>("filter_by_index1", 0, "specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line", false, "");
    cmd.add<string>("filter_by_index2", 0, "specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line", false, "");
    cmd.add<int>("filter_by_index_threshold", 0, "the allowed difference of index barcode for index filtering, default 0 means completely identical.", false, 0);
    
    // base correction in overlapped regions of paired end data
    cmd.add("correction", 'c', "enable base correction in overlapped regions (only for PE data), default is disabled");
    cmd.add<int>("overlap_len_require", 0, "the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default.", false, 30);
    cmd.add<int>("overlap_diff_limit", 0, "the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default.", false, 5);
    cmd.add<int>("overlap_diff_percent_limit", 0, "the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%.", false, 20);

    // umi
    cmd.add("umi", 'U', "enable unique molecular identifier (UMI) preprocessing");
    cmd.add<string>("umi_loc", 0, "specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none", false, "");
    cmd.add<int>("umi_len", 0, "if the UMI is in read1/read2, its length should be provided", false, 0);
    cmd.add<string>("umi_prefix", 0, "if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default", false, "");
    cmd.add<int>("umi_skip", 0, "if the UMI is in read1/read2, fastp can skip several bases following UMI, default is 0", false, 0);
    cmd.add<string>("umi_delim", 0, "delimiter to use between the read name and the UMI, default is :", false, ":");

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

    // deprecated options
    cmd.add("cut_by_quality5", 0, "DEPRECATED, use --cut_front instead.");
    cmd.add("cut_by_quality3", 0, "DEPRECATED, use --cut_tail instead.");
    cmd.add("cut_by_quality_aggressive", 0, "DEPRECATED, use --cut_right instead.");
    cmd.add("discard_unmerged", 0, "DEPRECATED, no effect now, see the introduction for merging.");
    
    cmd.parse_check(argc, argv);

    if(argc == 1) {
        cerr << cmd.usage() <<endl;
    }

    if(argc == 1) {
        //output citation information
        cerr << "Citation:" <<endl;
        cerr << "Chen, S. (2025). Fastp 1.0: An ultra‐fast all‐round tool for FASTQ data quality control and preprocessing. Imeta, 4(5), e70078." << endl;
        cerr << endl;
        return 0;
    }

    if(cmd.exist("discard_unmerged")) {
        cerr << "DEPRECATED: --discard_unmerged has no effect now, see the introduction for merging." << endl;
    }

    Options opt;

    // I/O
    opt.in1 = cmd.get<string>("in1");
    opt.in2 = cmd.get<string>("in2");
    opt.out1 = cmd.get<string>("out1");
    opt.out2 = cmd.get<string>("out2");
    opt.unpaired1 = cmd.get<string>("unpaired1");
    opt.unpaired2 = cmd.get<string>("unpaired2");
    opt.failedOut = cmd.get<string>("failed_out");
    opt.overlappedOut = cmd.get<string>("overlapped_out");
    // write to the same file
    if(opt.unpaired2.empty())
        opt.unpaired2 = opt.unpaired1;
    opt.compression = cmd.get<int>("compression");
    opt.readsToProcess = cmd.get<int>("reads_to_process");
    opt.phred64 = cmd.exist("phred64");
    opt.dontOverwrite = cmd.exist("dont_overwrite");
    opt.inputFromSTDIN = cmd.exist("stdin");
    opt.outputToSTDOUT = cmd.exist("stdout");
    opt.interleavedInput = cmd.exist("interleaved_in");
    opt.verbose = cmd.exist("verbose");
    opt.fixMGI = cmd.exist("fix_mgi_id");

    // duplication evaluation and deduplication
    opt.duplicate.dedup = cmd.exist("dedup");
    opt.duplicate.enabled = !cmd.exist("dont_eval_duplication") || cmd.exist("dedup") ;
    if(!cmd.exist("dup_calc_accuracy")) {
        if(opt.duplicate.dedup)
            opt.duplicate.accuracyLevel = 3;
        else
            opt.duplicate.accuracyLevel = 1;
    } else {
        opt.duplicate.accuracyLevel = min(6, max(1, cmd.get<int>("dup_calc_accuracy")));
    }

    // merge PE
    opt.merge.enabled = cmd.exist("merge");
    opt.merge.out = cmd.get<string>("merged_out");
    opt.merge.includeUnmerged = cmd.exist("include_unmerged");

    // adapter cutting
    opt.adapter.enabled = !cmd.exist("disable_adapter_trimming");
    opt.adapter.detectAdapterForPE = cmd.exist("detect_adapter_for_pe");
    opt.adapter.allowGapOverlapTrimming = cmd.exist("allow_gap_overlap_trimming");
    opt.adapter.dimerMaxLen = cmd.get<int>("dimer_max_len");
    opt.adapter.sequence = cmd.get<string>("adapter_sequence");
    opt.adapter.sequenceR2 = cmd.get<string>("adapter_sequence_r2");
    opt.adapter.fastaFile = cmd.get<string>("adapter_fasta");
    if(opt.adapter.sequenceR2=="auto" && !opt.adapter.detectAdapterForPE && opt.adapter.sequence != "auto") {
        opt.adapter.sequenceR2 = opt.adapter.sequence;
    }
    if(!opt.adapter.fastaFile.empty()) {
        opt.loadFastaAdapters();
    }

    // trimming
    opt.trim.front1 = cmd.get<int>("trim_front1");
    opt.trim.tail1 = cmd.get<int>("trim_tail1");
    opt.trim.maxLen1 = cmd.get<int>("max_len1");
    // read2 settings follows read1 if it's not specified
    if(cmd.exist("trim_front2"))
        opt.trim.front2 = cmd.get<int>("trim_front2");
    else
        opt.trim.front2 = opt.trim.front1;
    if(cmd.exist("trim_tail2"))
        opt.trim.tail2 = cmd.get<int>("trim_tail2");
    else
        opt.trim.tail2 = opt.trim.tail1;
    if(cmd.exist("max_len2"))
        opt.trim.maxLen2 = cmd.get<int>("max_len2");
    else
        opt.trim.maxLen2 = opt.trim.maxLen1;

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
    opt.qualityCut.enabledFront = cmd.exist("cut_front");
    // back compatible with old versions
    if(!opt.qualityCut.enabledFront){
        opt.qualityCut.enabledFront = cmd.exist("cut_by_quality5");
        if(opt.qualityCut.enabledFront)
            cerr << "WARNING: cut_by_quality5 is deprecated, please use cut_front instead." << endl;
    }
    opt.qualityCut.enabledTail = cmd.exist("cut_tail");
    // back compatible with old versions
    if(!opt.qualityCut.enabledFront){
        opt.qualityCut.enabledFront = cmd.exist("cut_by_quality3");
        if(opt.qualityCut.enabledFront)
            cerr << "WARNING: cut_by_quality3 is deprecated, please use cut_tail instead." << endl;
    }
    opt.qualityCut.enabledRight = cmd.exist("cut_right");
    // back compatible with old versions
    if(!opt.qualityCut.enabledRight){
        opt.qualityCut.enabledRight = cmd.exist("cut_by_quality_aggressive");
        if(opt.qualityCut.enabledRight)
            cerr << "WARNING: cut_by_quality_aggressive is deprecated, please use cut_right instead." << endl;
    }

    opt.qualityCut.windowSizeShared = cmd.get<int>("cut_window_size");
    opt.qualityCut.qualityShared = cmd.get<int>("cut_mean_quality");

    if(cmd.exist("cut_front_window_size"))
        opt.qualityCut.windowSizeFront = cmd.get<int>("cut_front_window_size");
    else
        opt.qualityCut.windowSizeFront = opt.qualityCut.windowSizeShared;
    if(cmd.exist("cut_front_mean_quality"))
        opt.qualityCut.qualityFront = cmd.get<int>("cut_front_mean_quality");
    else
        opt.qualityCut.qualityFront = opt.qualityCut.qualityShared;

    if(cmd.exist("cut_tail_window_size"))
        opt.qualityCut.windowSizeTail = cmd.get<int>("cut_tail_window_size");
    else
        opt.qualityCut.windowSizeTail = opt.qualityCut.windowSizeShared;
    if(cmd.exist("cut_tail_mean_quality"))
        opt.qualityCut.qualityTail = cmd.get<int>("cut_tail_mean_quality");
    else
        opt.qualityCut.qualityTail = opt.qualityCut.qualityShared;

    if(cmd.exist("cut_right_window_size"))
        opt.qualityCut.windowSizeRight = cmd.get<int>("cut_right_window_size");
    else
        opt.qualityCut.windowSizeRight = opt.qualityCut.windowSizeShared;
    if(cmd.exist("cut_right_mean_quality"))
        opt.qualityCut.qualityRight = cmd.get<int>("cut_right_mean_quality");
    else
        opt.qualityCut.qualityRight = opt.qualityCut.qualityShared;

    // raise a warning if cutting option is not enabled but -W/-M is enabled
    if(!opt.qualityCut.enabledFront && !opt.qualityCut.enabledTail && !opt.qualityCut.enabledRight) {
        if(cmd.exist("cut_window_size") || cmd.exist("cut_mean_quality") 
            || cmd.exist("cut_front_window_size") || cmd.exist("cut_front_mean_quality") 
            || cmd.exist("cut_tail_window_size") || cmd.exist("cut_tail_mean_quality") 
            || cmd.exist("cut_right_window_size") || cmd.exist("cut_right_mean_quality"))
            cerr << "WARNING: you specified the options for cutting by quality, but forgot to enable any of cut_front/cut_tail/cut_right. This will have no effect." << endl;
    }

    // quality filtering
    opt.qualfilter.enabled = !cmd.exist("disable_quality_filtering");
    opt.qualfilter.qualifiedQual = num2qual(cmd.get<int>("qualified_quality_phred"));
    opt.qualfilter.unqualifiedPercentLimit = cmd.get<int>("unqualified_percent_limit");
    opt.qualfilter.avgQualReq = cmd.get<int>("average_qual");
    opt.qualfilter.nBaseLimit = cmd.get<int>("n_base_limit");

    // length filtering
    opt.lengthFilter.enabled = !cmd.exist("disable_length_filtering");
    opt.lengthFilter.requiredLength = cmd.get<int>("length_required");
    opt.lengthFilter.maxLength = cmd.get<int>("length_limit");

    // low complexity filter
    opt.complexityFilter.enabled = cmd.exist("low_complexity_filter");
    opt.complexityFilter.threshold = (min(100, max(0, cmd.get<int>("complexity_threshold")))) / 100.0;

    // overlap correction
    opt.correction.enabled = cmd.exist("correction");
    opt.overlapRequire = cmd.get<int>("overlap_len_require");
    opt.overlapDiffLimit = cmd.get<int>("overlap_diff_limit");
    opt.overlapDiffPercentLimit = cmd.get<int>("overlap_diff_percent_limit");

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

    if(opt.inputFromSTDIN || opt.in1=="/dev/stdin") {
        if(opt.split.needEvaluation) {
            error_exit("Splitting by file number is not supported in STDIN mode");
        }
    }

    // umi
    opt.umi.enabled = cmd.exist("umi");
    opt.umi.length = cmd.get<int>("umi_len");
    opt.umi.prefix = cmd.get<string>("umi_prefix");
    opt.umi.skip = cmd.get<int>("umi_skip");
    opt.umi.delimiter = cmd.get<string>("umi_delim");
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

    // filtering by index
    string blacklist1 = cmd.get<string>("filter_by_index1");
    string blacklist2 = cmd.get<string>("filter_by_index2");
    int indexFilterThreshold = cmd.get<int>("filter_by_index_threshold");
    opt.initIndexFiltering(blacklist1, blacklist2, indexFilterThreshold);

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    time_t t1 = time(NULL);

    bool supportEvaluation = !opt.inputFromSTDIN && opt.in1!="/dev/stdin";

    Evaluator eva(&opt);
    if(supportEvaluation) {
        eva.evaluateSeqLen();

        if(opt.overRepAnalysis.enabled)
            eva.evaluateOverRepSeqs();
    }

    long readNum = 0;

    // using evaluator to guess how many reads in total
    if(opt.shallDetectAdapter(false)) {
        if(!supportEvaluation)
            cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
        else {
            cerr << "Detecting adapter sequence for read1..." << endl;
            string adapt = eva.evalAdapterAndReadNum(readNum, false);
            if(adapt.length() > 60 )
                adapt.resize(0, 60);
            if(adapt.length() > 0 ) {
                opt.adapter.sequence = adapt;
                opt.adapter.detectedAdapter1 = adapt;
            } else {
                cerr << "No adapter detected for read1" << endl;
                opt.adapter.sequence = "";
            }
            cerr << endl;
        }
    }
    if(opt.shallDetectAdapter(true)) {
        if(!supportEvaluation)
            cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
        else {
            cerr << "Detecting adapter sequence for read2..." << endl;
            string adapt = eva.evalAdapterAndReadNum(readNum, true);
            if(adapt.length() > 60 )
                adapt.resize(0, 60);
            if(adapt.length() > 0 ) {
                opt.adapter.sequenceR2 = adapt;
                opt.adapter.detectedAdapter2 = adapt;
            } else {
                cerr << "No adapter detected for read2" << endl;
                opt.adapter.sequenceR2 = "";
            }
            cerr << endl;
        }
    }

    opt.validate();

    // using evaluator to guess how many reads in total
    if(opt.split.needEvaluation && supportEvaluation) {
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
    if(!cmd.exist("trim_poly_g") && !cmd.exist("disable_trim_poly_g") && supportEvaluation) {
        bool twoColorSystem = eva.isTwoColorSystem();
        if(twoColorSystem){
            opt.polyGTrim.enabled = true;
        }
    }

    Processor p(&opt);
    p.process();
    
    time_t t2 = time(NULL);

    cerr << endl << "JSON report: " << opt.jsonFile << endl;
    cerr << "HTML report: " << opt.htmlFile << endl;
    cerr << endl << command << endl;
    cerr << "fastp v" << FASTP_VER << ", time used: " << (t2)-t1 << " seconds" << endl;

    return 0;
}
