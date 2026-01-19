#include "adaptertrimmer.h"
#include "matcher.h"

AdapterTrimmer::AdapterTrimmer(){
}


AdapterTrimmer::~AdapterTrimmer(){
}

bool AdapterTrimmer::trimByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr, int diffLimit, int overlapRequire, double diffPercentLimit) {
    OverlapResult ov = OverlapAnalysis::analyze(r1, r2, diffLimit, overlapRequire, diffPercentLimit);
    return trimByOverlapAnalysis(r1, r2, fr, ov);
}

bool AdapterTrimmer::trimByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr, OverlapResult ov, int frontTrimmed1, int frontTrimmed2) {
    int ol = ov.overlap_len;
    if(ov.overlapped && ov.offset < 0) {

        //5'      ......frontTrimmed1......|------------------------------------------|----- 3'
        //3' -----|-------------------------------------------|......frontTrimmed2.....      5'

        int len1 = min(r1->length(), ol + frontTrimmed2);
        int len2 = min(r2->length(), ol + frontTrimmed1);
        string adapter1 = r1->mSeq->substr(len1, r1->length() - len1);
        string adapter2 = r2->mSeq->substr(len2, r2->length() - len2);

        if(_DEBUG) {
            cerr << adapter1 << endl;
            cerr << adapter2 << endl;
            cerr << "frontTrimmed2: " << frontTrimmed1 << endl;
            cerr << "frontTrimmed2: " << frontTrimmed2 << endl;
            cerr << "overlap:" << ov.offset << "," << ov.overlap_len << ", " << ov.diff << endl;
            r1->print();
            r2->reverseComplement()->print();
            cerr <<endl;
        }
        r1->resize(len1);
        r2->resize(len2);

        fr->addAdapterTrimmed(adapter1, adapter2);
        return true;
    }
    return false;
}

bool AdapterTrimmer::trimByMultiSequences(Read* r, FilterResult* fr, vector<string>& adapterList, bool isR2, bool incTrimmedCounter) {
    int matchReq = 4;
    if(adapterList.size() > 16)
        matchReq = 5;
    if(adapterList.size() > 256)
        matchReq = 6;
    bool trimmed = false;

    string* originalSeq = r->mSeq;
    for(int i=0; i<adapterList.size(); i++) {
        trimmed |= trimBySequence(r, NULL, adapterList[i], isR2, matchReq);
    }

    if(trimmed) {
        string adapter = originalSeq->substr(r->length(), originalSeq->length() - r->length());
        if(fr)
            fr->addAdapterTrimmed(adapter, isR2, incTrimmedCounter);
        else
            cerr << adapter << endl;
    }

    return trimmed;
}

bool AdapterTrimmer::trimBySequence(Read* r, FilterResult* fr, string& adapterseq, bool isR2, int matchReq) {
    const int allowOneMismatchForEach = 8;

    int rlen = r->length();
    int alen = adapterseq.length();

    const char* adata = adapterseq.c_str();
    const char* rdata = r->mSeq->c_str();

    if(alen < matchReq)
        return false;

    int pos=0;
    bool found = false;
    int start = 0;
    if(alen >= 16)
        start = -4;
    else if(alen >= 12)
        start = -3;
    else if(alen >= 8)
        start = -2;
    // we start from negative numbers since the Illumina adapter dimer usually have the first A skipped as A-tailing
    // try exact match with hamming distance (no insertion of deletion)
    for(pos = start; pos<rlen-matchReq; pos++) {
        int cmplen = min(rlen - pos, alen);
        int allowedMismatch = cmplen/allowOneMismatchForEach;
        int mismatch = 0;
        bool matched = true;
        for(int i=max(0, -pos); i<cmplen; i++) {
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

    // if failed to exact match, we try one gap
    // to lower computational cost, we only allow one gap, and it's much enough for short reads
    // we try insertion in the sequence
    bool hasInsertion = false;
    if(!found) {
        for(pos = 0; pos<rlen-matchReq-1; pos++) {
            int cmplen = min(rlen - pos - 1, alen);
            int allowedMismatch = cmplen/allowOneMismatchForEach -1;
            bool matched = Matcher::matchWithOneInsertion(rdata, adata, cmplen, allowedMismatch);
            if(matched) {
                found = true;
                hasInsertion = true;
                //cerr << ".";
                break;
            }
        }
    }

    // if failed to exact match, and failed to match with one insertion in sequence
    // we then try deletion in the sequence
    bool hasDeletion = false;
    if(!found) {
        for(pos = 0; pos<rlen-matchReq; pos++) {
            int cmplen = min(rlen - pos, alen - 1);
            int allowedMismatch = cmplen/allowOneMismatchForEach -1;
            bool matched = Matcher::matchWithOneInsertion(adata, rdata, cmplen, allowedMismatch);
            if(matched) {
                found = true;
                hasDeletion = true;
                //cerr << "|";
                break;
            }
        }
    }

    if(found) {
        if(pos < 0) {
            string adapter = adapterseq.substr(0, alen+pos);
            r->mSeq->resize(0);
            r->mQuality->resize(0);
            if(fr) {
                fr->addAdapterTrimmed(adapter, isR2);
            }

        } else {
            string adapter = r->mSeq->substr(pos, rlen-pos);
            r->resize(pos);
            if(fr) {
                fr->addAdapterTrimmed(adapter, isR2);
            }
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
    if (*r.mSeq != "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA")
        return false;

    Read read("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGGAAATTTCCCGGGAAATTTCCCGGGATCGATCGATCGATCGAATTCC",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");
    vector<string> adapterList;
    adapterList.push_back("GCTAGCTAGCTAGCTA");
    adapterList.push_back("AAATTTCCCGGGAAATTTCCCGGG");
    adapterList.push_back("ATCGATCGATCGATCG");
    adapterList.push_back("AATTCCGGAATTCCGG");
    trimmed = AdapterTrimmer::trimByMultiSequences(&read, NULL, adapterList);
    if (*read.mSeq != "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG") {
        cerr << read.mSeq << endl;
        return false;
    }

    return true;
}