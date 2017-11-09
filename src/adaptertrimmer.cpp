#include "adaptertrimmer.h"

AdapterTrimmer::AdapterTrimmer(){
}


AdapterTrimmer::~AdapterTrimmer(){
}

bool AdapterTrimmer::trimByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr) {
    OverlapResult ov = OverlapAnalysis::analyze(r1, r2);
    int ol = ov.overlap_len;
    if(ov.diff<=5 && ov.overlapped && ov.offset < 0 && ol > r1->length()/2) {
        string adapter1 = r1->mSeq.mStr.substr(ol, r1->length() - ol);
        string adapter2 = r2->mSeq.mStr.substr(ol, r2->length() - ol);

        if(_DEBUG) {
            cout << adapter1 << endl;
            cout << adapter2 << endl;
            cout << "overlap:" << ov.offset << "," << ov.overlap_len << ", " << ov.diff << endl;
            r1->print();
            r2->reverseComplement()->print();
            cout <<endl;
        }

        r1->mSeq.mStr = r1->mSeq.mStr.substr(0, ol);
        r1->mQuality = r1->mQuality.substr(0, ol);
        r2->mSeq.mStr = r2->mSeq.mStr.substr(0, ol);
        r2->mQuality = r2->mQuality.substr(0, ol);

        fr->addAdapterTrimmed(adapter1, adapter2);
        return true;
    }
    return false;
}

bool AdapterTrimmer::trimBySequence(Read* r1, FilterResult* fr, string& adapterseq) {
    const int matchReq = 4;
    const int allowOneMismatchForEach = 8;

    int rlen = r1->length();
    int alen = adapterseq.length();

    const char* adata = adapterseq.c_str();
    const char* rdata = r1->mSeq.mStr.c_str();

    if(alen < matchReq)
        return false;

    int l=0;
    bool found = false;
    for(l = matchReq; l<min(rlen, alen); l++) {
        int allowedMismatch = l/allowOneMismatchForEach;
        int mismatch = 0;
        bool matched = true;
        for(int i=0; i<l; i++) {
            if( adata[i] != rdata[rlen - l + i] ){
                mismatch++;
                if(mismatch > allowedMismatch) {
                    matched = false;
                    break;
                }
            }
        }
        if(matched) {
            found = true;
            break;
        }

    }

    if(found) {
        string adapter1 = r1->mSeq.mStr.substr(rlen - l, l);
        r1->mSeq.mStr = r1->mSeq.mStr.substr(0, rlen - l);
        r1->mQuality = r1->mQuality.substr(0, rlen - l);
        if(fr) {
            fr->addAdapterTrimmed(adapter1);
        }
        return true;
    }

    return false;
}

bool AdapterTrimmer::test() {
    Read r("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////E");
    string adapter = "TTTTCCACGGGGATACTACTG";
    bool trimmed = AdapterTrimmer::trimBySequence(&r, NULL, adapter);
    return r.mSeq.mStr == "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA";
}