#include <stdio.h>
#include "fastqreader.h"
#include "unittest.h"
#include <time.h>
#include "cmdline.h"
#include <sstream>
#include "util.h"
#include "options.h"
#include "processor.h"

string command;

int main(int argc, char* argv[]){
    if (argc == 2 && strcmp(argv[1], "test")==0){
        UnitTest tester;
        tester.run();
        return 0;
    }
    cmdline::parser cmd;
    // input/output
    cmd.add<string>("in1", 'i', "read1 input file name", true, "");
    cmd.add<string>("out1", 'o', "read1 output file name", false, "");
    cmd.add<string>("in2", 'I', "read2 input file name", false, "");
    cmd.add<string>("out2", 'O', "read2 output file name", false, "");
    cmd.add<int>("compression", 'z', "compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 2.", false, 2);

    // trimming
    cmd.add<int>("trim_front1", 'f', "trimming how many bases in front for read1, default is 0", false, 0);
    cmd.add<int>("trim_tail1", 't', "trimming how many bases in tail for read1, default is 0", false, 0);
    cmd.add<int>("trim_front2", 'F', "trimming how many bases in front for read2. If it's not specified, it will follow read1's settings", false, 0);
    cmd.add<int>("trim_tail2", 'T', "trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings", false, 0);

    // quality filtering
    cmd.add("disable_quality_filtering", 'Q', "quality filtering is enabled by default. If this option is enabled, quality filtering is disabled");
    cmd.add<int>("qualified_quality_phred", 'q', "the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified.", false, 15);
    cmd.add<int>("unqualified_percent_limit", 'u', "how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%", false, 40);
    cmd.add<int>("n_base_limit", 'n', "if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5", false, 5);

    // length filtering
    cmd.add<int>("length_required", 'l', "length filtering will be enabled if this argument is specified, reads shorter than length_required will be discarded.", false, 30);
    

    // reporting
    cmd.add<string>("json", 'j', "the json format report file name", false, "fastp.json");
    cmd.add<string>("html", 'h', "the html format report file name", false, "fastp.html");

    // threading
    cmd.add<int>("thread", 'w', "worker thread number, default is 3", false, 3);

    cmd.parse_check(argc, argv);

    Options opt;

    // I/O
    opt.in1 = cmd.get<string>("in1");
    opt.in2 = cmd.get<string>("in2");
    opt.out1 = cmd.get<string>("out1");
    opt.out2 = cmd.get<string>("out2");
    opt.compression = cmd.get<int>("compression");

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

    // quality filtering
    opt.qualfilter.enabled = !cmd.exist("disable_quality_filtering");
    opt.qualfilter.qualifiedQual = num2qual(cmd.get<int>("qualified_quality_phred"));
    opt.qualfilter.unqualifiedPercentLimit = cmd.get<int>("unqualified_percent_limit");
    opt.qualfilter.nBaseLimit = cmd.get<int>("n_base_limit");

    // length filtering
    opt.lengthFilter.enabled = cmd.exist("length_required");
    opt.lengthFilter.requiredLength = cmd.get<int>("length_required");

    // threading
    opt.thread = cmd.get<int>("thread");

    // reporting
    opt.jsonFile = cmd.get<string>("json");
    opt.htmlFile = cmd.get<string>("html");

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    opt.validate();

    time_t t1 = time(NULL);

    Processor p(&opt);
    p.process();
    
    time_t t2 = time(NULL);

    cout << endl << "JSON report: " << opt.jsonFile << endl;
    cout << "HTML report: " << opt.htmlFile << endl;
    cout << endl << command << endl;
    cout << "fastp v" << FASTP_VER << ", time used: " << (t2)-t1 << " seconds" << endl;
}