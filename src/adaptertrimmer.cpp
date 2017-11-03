#include "adaptertrimmer.h"

AdapterTrimmer::AdapterTrimmer(){
}


AdapterTrimmer::~AdapterTrimmer(){
}

bool AdapterTrimmer::trimByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr) {
    OverlapResult ov = OverlapAnalysis::analyze(r1, r2);
    int ol = ov.overlap_len;
    if(ov.overlapped && ov.offset < 0) {
        string adapter1 = r1->mSeq.mStr.substr(ol, r1->length() - ol);
        string adapter2 = r2->mSeq.mStr.substr(ol, r2->length() - ol);

        r1->mSeq.mStr = r1->mSeq.mStr.substr(0, ol);
        r1->mQuality = r1->mQuality.substr(0, ol);
        r2->mSeq.mStr = r2->mSeq.mStr.substr(0, ol);
        r2->mQuality = r2->mQuality.substr(0, ol);

        fr->addAdapterTrimmed(adapter1, adapter2);
        return true;
    }
    return false;
}