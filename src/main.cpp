#include <stdio.h>
#include "fastqreader.h"
#include "unittest.h"
#include <time.h>
#include "cmdline.h"
#include <sstream>
#include "util.h"
#include "options.h"

string command;

int main(int argc, char* argv[]){
    if (argc == 2 && strcmp(argv[1], "test")==0){
        UnitTest tester;
        tester.run();
        return 0;
    }
    cmdline::parser cmd;
    cmd.add<string>("in1", 'i', "read1 input file name", true, "");
    cmd.add<string>("out1", 'o', "read1 output file name", true, "");
    cmd.add<string>("in2", 'I', "read2 input file name", false, "");
    cmd.add<string>("out2", 'O', "read2 output file name", false, "");
    cmd.add<string>("json", 'j', "the json format report file name", false, "fastp.json");
    cmd.add<string>("html", 'h', "the html format report file name", false, "fastp.html");
    cmd.add<int>("thread", 't', "worker thread number, default is 4", false, 4);
    cmd.parse_check(argc, argv);

    Options opt;

    opt.in1 = cmd.get<string>("in1");
    opt.in2 = cmd.get<string>("in2");
    opt.out1 = cmd.get<string>("out1");
    opt.out2 = cmd.get<string>("out2");
    opt.thread = cmd.get<int>("thread");

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    check_file_valid(opt.in1);
    if(opt.in2 != "")
        check_file_valid(opt.in2);

    time_t t1 = time(NULL);
    
    
    time_t t2 = time(NULL);
    printf("\n%s\n", command.c_str());
    printf("fastp v%s, time used: %ld seconds\n", FASTP_VER, (t2-t1));
}