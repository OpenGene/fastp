#include "overlapanalysis.h"

OverlapAnalysis::OverlapAnalysis(){
}


OverlapAnalysis::~OverlapAnalysis(){
}

// ported from the python code of AfterQC
OverlapResult OverlapAnalysis::analyze(Sequence& r1, Sequence& r2) {
    Sequence rcr2 = ~r2;
    int len1 = r1.length();
    int len2 = rcr2.length();
    // use the pointer directly for speed
    const char* str1 = r1.mStr.c_str();
    const char* str2 = rcr2.mStr.c_str();

    int limit_distance = 3;
    int overlap_require = 30;
    int complete_compare_require = 50;

    int overlap_len = 0;
    int offset = 0;
    int diff = 0;

    // forward
    // a match of less than 30 is considered as unconfident
    while (offset < len1-overlap_require) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1 - offset, len2);

        diff = 0;
        int i = 0;
        for (i=0; i<overlap_len; i++) {
            if (str1[offset + i] != str2[i]){
                diff += 1;
                if (diff >= limit_distance && i < complete_compare_require)
                    break;
            }
        }
        
        if (diff < limit_distance || (diff >= limit_distance && i>complete_compare_require)){
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
    while (offset > -(len2-overlap_require)){
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1,  len2- abs(offset));

        diff = 0;
        int i = 0;
        for (i=0; i<overlap_len; i++) {
            if (str1[i] != str2[-offset + i]){
                diff += 1;
                if (diff >= limit_distance && i < complete_compare_require)
                    break;
            }
        }
        
        if (diff < limit_distance || (diff >= limit_distance && i>complete_compare_require)){
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

bool OverlapAnalysis::test(){
    //Sequence r1("CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGCCGCTGGAGGTCTCCC");
    //Sequence r2("ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGCCCGTAGGCGCGGCTCCC");

    Sequence r1("CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGC");
    Sequence r2("ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGTCC");
    
    OverlapResult ov = OverlapAnalysis::analyze(r1, r2);

    return ov.overlapped && ov.offset == 10 && ov.overlap_len == 79 && ov.diff == 1;
}