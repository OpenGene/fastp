#include "options.h"
#include "util.h"

Options::Options(){
    in1 = "";
    in2 = "";
    out1 = "";
    out2 = "";
    thread = 1;
    compression = 2;
}

void Options::init() {
}

bool Options::isPaired() {
    return in2.length() > 0;
}

bool Options::adapterCuttingEnabled() {
    return adapter.enabled && isPaired();
}

bool Options::validate() {
    if(in1.empty()) {
        error_exit("read1 input should be specified by --in1");
    } else {
        check_file_valid(in1);
    }

    if(in2.empty() && !out2.empty()) {
        error_exit("read2 output is specified (--out2), but read2 input is not specified (--in2)");
    }

    if(!in2.empty()) {
        check_file_valid(in2);

        if(!out1.empty() && out2.empty()) {
            error_exit("paired-end input, read1 output should be specified together with read2 output (--out2 needed) ");
        }
        if(out1.empty() && !out2.empty()) {
            error_exit("paired-end input, read1 output should be specified (--out1 needed) together with read2 output ");
        }
    }

    if(compression < 1 || compression > 9)
        error_exit("compression level (--compression) should be between 1 ~ 9, 1 for fastest, 9 for smallest");

    if(thread < 1 || thread > 16)
        error_exit("thread number (--thread) should be 1 ~ 16, suggest 1 ~ 8");

    if(trim.front1 < 0 || trim.front1 > 30)
        error_exit("trim_front1 (--trim_front1) should be 0 ~ 30, suggest 0 ~ 4");

    if(trim.tail1 < 0 || trim.tail1 > 30)
        error_exit("trim_tail1 (--trim_tail1) should be 0 ~ 30, suggest 0 ~ 4");

    if(trim.front2 < 0 || trim.front2 > 30)
        error_exit("trim_front2 (--trim_front2) should be 0 ~ 30, suggest 0 ~ 4");

    if(trim.tail2 < 0 || trim.tail2 > 30)
        error_exit("trim_tail2 (--trim_tail2) should be 0 ~ 30, suggest 0 ~ 4");

    if(qualfilter.qualifiedQual - 33 < 0 || qualfilter.qualifiedQual - 33 > 50)
        error_exit("qualitified phred (--qualified_quality_phred) should be 0 ~ 50, suggest 10 ~ 20");

    if(qualfilter.unqualifiedPercentLimit < 0 || qualfilter.unqualifiedPercentLimit > 100)
        error_exit("unqualified percent limit (--unqualified_percent_limit) should be 0 ~ 100, suggest 20 ~ 60");

    if(qualfilter.nBaseLimit < 0 || qualfilter.nBaseLimit > 50)
        error_exit("N base limit (--n_base_limit) should be 0 ~ 50, suggest 3 ~ 10");

    if(lengthFilter.requiredLength < 0 )
        error_exit("length requirement (--length_required) should be >0, suggest 20 ~ 100");

    if(split.enabled ) {
        if(split.number < 2 || split.number > 1000)
            error_exit("you have enabled splitting output to multiple files, the number of files (--split) should be 2 ~ 1000.");
        if(split.digits < 0 || split.digits > 10)
            error_exit("you have enabled splitting output to multiple files, the digits number of file name prefix (--split_prefix_digits) should be 0 ~ 10.");
        if(split.number % thread != 0) {
            cerr << "you have enabled splitting output to multiple files, but the number of files (--split) should be a multiple of threads (--thread)." << endl;
            cerr << "current thread number setting is: " << thread << ", you can change the thread number by (--thread)" << endl;
            cerr << "or set the number of files (--split) to be a multiple of " << thread << ", for example: " << thread << ", " << thread*2 << ", " << thread*3 << "... " << endl;
            cerr << endl;
            exit(-1);
        }
    }

    return true;
}