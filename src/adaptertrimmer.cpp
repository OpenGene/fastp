#include "adaptertrimmer.h"

AdapterTrimmer::AdapterTrimmer(){
}


AdapterTrimmer::~AdapterTrimmer(){
}

bool AdapterTrimmer::trimByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr) {
    OverlapResult ov = OverlapAnalysis::analyze(r1, r2);
    return trimByOverlapAnalysis(r1, r2, fr, ov);
}

bool AdapterTrimmer::trimByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr, OverlapResult ov) {
    int ol = ov.overlap_len;
    if(ov.diff<=5 && ov.overlapped && ov.offset < 0 && ol > r1->length()/3) {
        string adapter1 = r1->mSeq.mStr.substr(ol, r1->length() - ol);
        string adapter2 = r2->mSeq.mStr.substr(ol, r2->length() - ol);

        if(_DEBUG) {
            cerr << adapter1 << endl;
            cerr << adapter2 << endl;
            cerr << "overlap:" << ov.offset << "," << ov.overlap_len << ", " << ov.diff << endl;
            r1->print();
            r2->reverseComplement()->print();
            cerr <<endl;
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

bool AdapterTrimmer::trimBySequence(Read* r, FilterResult* fr, string& adapterseq, int allowOneMismatchForEach, bool isR2) {
    const int matchReq = 4;

    int rlen = r->length();
    int alen = adapterseq.length();

    const char* adata = adapterseq.c_str();
    const char* rdata = r->mSeq.mStr.c_str();

    if(alen < matchReq)
        return false;

    int pos=0;
    bool found = false;
    for(pos = 0; pos<rlen-matchReq; pos++) {
        int cmplen = min(rlen - pos, alen);
        int allowedMismatch = cmplen/allowOneMismatchForEach;
        int mismatch = 0;
        bool matched = true;
        for(int i=0; i<cmplen; i++) {
            if( adata[i] != rdata[i+pos] ){
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
        string adapter = r->mSeq.mStr.substr(pos, rlen-pos);
        r->mSeq.mStr = r->mSeq.mStr.substr(0, pos);
        r->mQuality = r->mQuality.substr(0, pos);
        if(fr) {
            fr->addAdapterTrimmed(adapter, isR2);
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
    bool trimmed = AdapterTrimmer::trimBySequence(&r, NULL, adapter, 8);
    return r.mSeq.mStr == "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA";
}
