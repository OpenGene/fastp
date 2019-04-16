#include "overlapanalysis.h"

OverlapAnalysis::OverlapAnalysis(){
}


OverlapAnalysis::~OverlapAnalysis(){
}

OverlapResult OverlapAnalysis::analyze(Read* r1, Read* r2, int overlapDiffLimit, int overlapRequire, double diffPercentLimit) {
    return analyze(r1->mSeq, r2->mSeq, overlapDiffLimit, overlapRequire, diffPercentLimit);
}

// ported from the python code of AfterQC
OverlapResult OverlapAnalysis::analyze(Sequence& r1, Sequence& r2, int diffLimit, int overlapRequire, double diffPercentLimit) {
    Sequence rcr2 = ~r2;
    int len1 = r1.length();
    int len2 = rcr2.length();
    // use the pointer directly for speed
    const char* str1 = r1.mStr.c_str();
    const char* str2 = rcr2.mStr.c_str();

    int complete_compare_require = 50;

    int overlap_len = 0;
    int offset = 0;
    int diff = 0;

    // forward
    // a match of less than overlapRequire is considered as unconfident
    while (offset < len1-overlapRequire) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1 - offset, len2);
        int overlapDiffLimit = min(diffLimit, (int)(overlap_len * diffPercentLimit));

        diff = 0;
        int i = 0;
        for (i=0; i<overlap_len; i++) {
            if (str1[offset + i] != str2[i]){
                diff += 1;
                if (diff > overlapDiffLimit && i < complete_compare_require)
                    break;
            }
        }
        
        if (diff <= overlapDiffLimit || (diff > overlapDiffLimit && i>complete_compare_require)){
            OverlapResult ov;
            ov.overlapped = true;
            ov.offset = offset;
            ov.overlap_len = overlap_len;
            ov.diff = diff;
            return ov;
        }

        offset += 1;
    }


    // reverse
    // in this case, the adapter is sequenced since TEMPLATE_LEN < SEQ_LEN
    // check if distance can get smaller if offset goes negative
    // this only happens when insert DNA is shorter than sequencing read length, and some adapter/primer is sequenced but not trimmed cleanly
    // we go reversely
    offset = 0;
    while (offset > -(len2-overlapRequire)){
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1,  len2- abs(offset));
        int overlapDiffLimit = min(diffLimit, (int)(overlap_len * diffPercentLimit));

        diff = 0;
        int i = 0;
        for (i=0; i<overlap_len; i++) {
            if (str1[i] != str2[-offset + i]){
                diff += 1;
                if (diff > overlapDiffLimit && i < complete_compare_require)
                    break;
            }
        }
        
        if (diff <= overlapDiffLimit || (diff > overlapDiffLimit && i>complete_compare_require)){
            OverlapResult ov;
            ov.overlapped = true;
            ov.offset = offset;
            ov.overlap_len = overlap_len;
            ov.diff = diff;
            return ov;
        }

        offset -= 1;
    }

    OverlapResult ov;
    ov.overlapped = false;
    ov.offset = ov.overlap_len = ov.diff = 0;
    return ov;
}

Read* OverlapAnalysis::merge(Read* r1, Read* r2, OverlapResult ov) {
    int ol = ov.overlap_len;
    if(!ov.overlapped)
        return NULL;

    int len1 = ol + max(0, ov.offset);
    int len2 = 0; 
    if(ov.offset > 0)
        len2 = r2->length() - ol;

    Read* rr2 = r2->reverseComplement();
    string mergedSeq = r1->mSeq.mStr.substr(0, len1);
    if(ov.offset > 0) {
        mergedSeq += rr2->mSeq.mStr.substr(ol, len2);
    }

    string mergedQual = r1->mQuality.substr(0, len1);
    if(ov.offset > 0) {
        mergedQual += rr2->mQuality.substr(ol, len2);
    }

    delete rr2;

    string name = r1->mName + " merged_" + to_string(len1) + "_" + to_string(len2);
    Read* mergedRead = new Read(name, mergedSeq, r1->mStrand, mergedQual);

    return mergedRead;
}

bool OverlapAnalysis::test(){
    //Sequence r1("CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGCCGCTGGAGGTCTCCC");
    //Sequence r2("ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGCCCGTAGGCGCGGCTCCC");

    Sequence r1("CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGC");
    Sequence r2("ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGTCC");
    string qual1("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
    string qual2("#########################################################################################");
    
    OverlapResult ov = OverlapAnalysis::analyze(r1, r2, 2, 30, 0.2);

    Read read1("name1", r1, "+", qual1);
    Read read2("name2", r2, "+", qual2);

    Read* mergedRead = OverlapAnalysis::merge(&read1, &read2, ov);
    mergedRead->print();

    return ov.overlapped && ov.offset == 10 && ov.overlap_len == 79 && ov.diff == 1;
}