#include "jsonreporter.h"

JsonReporter::JsonReporter(Options* opt){
    mOptions = opt;
}


JsonReporter::~JsonReporter(){
}

void JsonReporter::report(FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2) {
    ofstream ofs;
    ofs.open(mOptions->jsonFile, ifstream::out);
    ofs << "{" << endl;

    if(result) {
        ofs << "\t" << "\"filtering_result\": " ;
        result -> reportJson(ofs, "\t");
    }

    if(preStats1) {
        ofs << "\t" << "\"read1_before_filtering\": " ;
        preStats1 -> reportJson(ofs, "\t");
    }

    if(postStats1) {
        ofs << "\t" << "\"read1_after_filtering\": " ;
        postStats1 -> reportJson(ofs, "\t");
    }

    if(preStats2) {
        ofs << "\t" << "\"read2_before_filtering\": " ;
        preStats2 -> reportJson(ofs, "\t");
    }

    if(postStats2) {
        ofs << "\t" << "\"read2_after_filtering\": " ;
        postStats2 -> reportJson(ofs, "\t");
    }

    ofs << "}";
}