#include "polyx.h"
#include "common.h"

PolyX::PolyX(){
}


PolyX::~PolyX(){
}

void PolyX::trimPolyG(Read* r1, Read* r2, FilterResult* fr, int compareReq, int allowOneMismatchForEach, int maxMismatch) {
    trimPolyG(r1, fr, compareReq, allowOneMismatchForEach, maxMismatch);
    trimPolyG(r2, fr, compareReq, allowOneMismatchForEach, maxMismatch);
}

void PolyX::trimPolyG(Read* r, FilterResult* fr, int compareReq, int allowOneMismatchForEach, int maxMismatch) {
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

void PolyX::trimPolyX(Read* r1, Read* r2, FilterResult* fr, int compareReq, int allowOneMismatchForEach, int maxMismatch) {
    trimPolyX(r1, fr, compareReq, allowOneMismatchForEach, maxMismatch);
    trimPolyX(r2, fr, compareReq, allowOneMismatchForEach, maxMismatch);
}

void PolyX::trimPolyX(Read* r, FilterResult* fr, int compareReq, int allowOneMismatchForEach, int maxMismatch) {
    const char* data = r->mSeq.mStr.c_str();

    int rlen = r->length();


    int atcgNumbers[4] = {0, 0, 0, 0};
    int pos = 0;
    for(pos=0; pos<rlen; pos++) {
        switch(data[rlen - pos - 1]) {
            case 'A':
                atcgNumbers[0]++;
                break;
            case 'T':
                atcgNumbers[1]++;
                break;
            case 'C':
                atcgNumbers[2]++;
                break;
            case 'G':
                atcgNumbers[3]++;
                break;
            case 'N':
                atcgNumbers[0]++;
                atcgNumbers[1]++;
                atcgNumbers[2]++;
                atcgNumbers[3]++;
                break;
            default:
                break;
        }

        int cmp = (pos+1);
        int allowedMismatch = min(maxMismatch, cmp/allowOneMismatchForEach);

        bool needToBreak = true;
        for(int b=0; b<4; b++) {
            if(cmp - atcgNumbers[b] <= allowedMismatch)
                needToBreak = false;
        }
        if(needToBreak && (pos >= allowOneMismatchForEach || pos+1 >= compareReq-1)) {
            break;
        }
    }

    // has polyX
    if(pos+1 >= compareReq) {
        // find the poly
        int poly;
        int maxCount = -1;
        for(int b=0; b<4; b++) {
            if(atcgNumbers[b] > maxCount){
                maxCount = atcgNumbers[b];
                poly = b;
            }
        }
        char polyBase = ATCG_BASES[poly];
        while(data[rlen - pos - 1] != polyBase && pos>=0)
            pos--;

        r->resize(rlen - pos - 1);
        if(fr)
          fr->addPolyXTrimmed(poly, pos + 1);
    }
}

bool PolyX::test() {

    Read r("@name",
        "ATTTTAAAAAAAAAATAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAT",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////E");

    FilterResult fr(NULL, false);
    PolyX::trimPolyX(&r, NULL, 10, 8, 5);
    r.print();

    return r.mSeq.mStr == "ATTTT" && fr.getTotalPolyXTrimmedReads() == 1 && fr.getTotalPolyXTrimmedBases() == 51;
}
