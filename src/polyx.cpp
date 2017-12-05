#include "polyx.h"

PolyX::PolyX(){
}


PolyX::~PolyX(){
}

void PolyX::trimPolyG(Read* r1, Read* r2, FilterResult* fr) {
    trimPolyG(r1, fr);
    trimPolyG(r2, fr);
}

void PolyX::trimPolyG(Read* r, FilterResult* fr) {
    const int compareReq = 10;
    const int allowOneMismatchForEach = 8;
    const int maxMismatch = 5;

    const char* data = r->mSeq.mStr.c_str();

    int rlen = r->length();

    int mismatch = 0;
    int i = 0;
    int firstGPos = rlen - 1;
    for(i=0; i< rlen; i++) {
        if(data[rlen - i - 1] != 'G') {
            mismatch++;
        } else {
            firstGPos = rlen - i -1;
        }

        int allowedMismatch = (i+1)/allowOneMismatchForEach;
        if(mismatch > maxMismatch || (mismatch>allowedMismatch && i>= compareReq-1) )
            break;
    }

    if(i >= compareReq) {
        r->resize(firstGPos);
    }
}

bool PolyX::test() {
    Read r("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAAGGGGGGGGGGGGGGGGGAGG",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////E");
    PolyX::trimPolyG(&r, NULL);
    r.print();
    return r.mSeq.mStr == "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAA";
}