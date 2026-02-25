#include "matcher.h"
#include "overlapanalysis.h"
#include "simd.h"

OverlapAnalysis::OverlapAnalysis(){
}


OverlapAnalysis::~OverlapAnalysis(){
}

OverlapResult OverlapAnalysis::analyze(Read* r1, Read* r2, int overlapDiffLimit, int overlapRequire, double diffPercentLimit, bool allowGap) {
    return analyze(r1->mSeq, r2->mSeq, overlapDiffLimit, overlapRequire, diffPercentLimit, allowGap);
}

// ported from the python code of AfterQC
OverlapResult OverlapAnalysis::analyze(string*  r1, string*  r2, int diffLimit, int overlapRequire, double diffPercentLimit, bool allowGap) {
    string rcr2 = Sequence::reverseComplement(r2);
    int len1 = r1->length();
    int len2 = rcr2.length();
    // use the pointer directly for speed
    const char* str1 = r1->c_str();
    const char* str2 = rcr2.c_str();

    int complete_compare_require = 50;

    int overlap_len = 0;
    int offset = 0;
    int diff = 0;

    // forward with no gap
    // a match of less than overlapRequire is considered as unconfident
    while (offset < len1-overlapRequire) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1 - offset, len2);
        int overlapDiffLimit = min(diffLimit, (int)(overlap_len * diffPercentLimit));

        diff = fastp_simd::countMismatches(str1 + offset, str2, overlap_len);

        // Accept if within limit, or for long overlaps where the prefix matches well
        bool accepted = (diff <= overlapDiffLimit);
        if (!accepted && overlap_len > complete_compare_require) {
            int prefixDiff = fastp_simd::countMismatches(str1 + offset, str2, complete_compare_require);
            accepted = (prefixDiff <= overlapDiffLimit);
        }

        if (accepted) {
            OverlapResult ov;
            ov.overlapped = true;
            ov.offset = offset;
            ov.overlap_len = overlap_len;
            ov.diff = diff;
            ov.hasGap = false;
            return ov;
        }

        offset += 1;
    }


    // reverse with no gap
    // in this case, the adapter is sequenced since TEMPLATE_LEN < SEQ_LEN
    // check if distance can get smaller if offset goes negative
    // this only happens when insert DNA is shorter than sequencing read length, and some adapter/primer is sequenced but not trimmed cleanly
    // we go reversely
    offset = 0;
    while (offset > -(len2-overlapRequire)){
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1,  len2- abs(offset));
        int overlapDiffLimit = min(diffLimit, (int)(overlap_len * diffPercentLimit));

        diff = fastp_simd::countMismatches(str1, str2 + (-offset), overlap_len);

        bool accepted = (diff <= overlapDiffLimit);
        if (!accepted && overlap_len > complete_compare_require) {
            int prefixDiff = fastp_simd::countMismatches(str1, str2 + (-offset), complete_compare_require);
            accepted = (prefixDiff <= overlapDiffLimit);
        }

        if (accepted) {
            OverlapResult ov;
            ov.overlapped = true;
            ov.offset = offset;
            ov.overlap_len = overlap_len;
            ov.diff = diff;
            ov.hasGap = false;
            return ov;
        }

        offset -= 1;
    }

    if(allowGap) {
        // forward with one gap
        offset = 0;
        while (offset < len1-overlapRequire) {
            // the overlap length of r1 & r2 when r2 is move right for offset
            overlap_len = min(len1 - offset, len2);
            int overlapDiffLimit = min(diffLimit, (int)(overlap_len * diffPercentLimit));

            int diff = Matcher::diffWithOneInsertion(str1 + offset, str2, overlap_len-1, overlapDiffLimit);
            if(diff <0 || diff > overlapDiffLimit)
                diff = Matcher::diffWithOneInsertion(str2, str1 + offset, overlap_len-1, overlapDiffLimit);
            
            if (diff <= overlapDiffLimit && diff >=0){
                OverlapResult ov;
                ov.overlapped = true;
                ov.offset = offset;
                ov.overlap_len = overlap_len;
                ov.diff = diff;
                ov.hasGap = true;
                return ov;
            }

            offset += 1;
        }

        // reverse with one gap
        offset = 0;
        while (offset > -(len2-overlapRequire)){
            // the overlap length of r1 & r2 when r2 is move right for offset
            overlap_len = min(len1,  len2- abs(offset));
            int overlapDiffLimit = min(diffLimit, (int)(overlap_len * diffPercentLimit));

            int diff = Matcher::diffWithOneInsertion(str1, str2-offset, overlap_len-1, overlapDiffLimit);
            if(diff <0 || diff > overlapDiffLimit)
                diff = Matcher::diffWithOneInsertion(str2-offset, str1, overlap_len-1, overlapDiffLimit);
            
            if (diff <= overlapDiffLimit && diff >=0){
                OverlapResult ov;
                ov.overlapped = true;
                ov.offset = offset;
                ov.overlap_len = overlap_len;
                ov.diff = diff;
                ov.hasGap = true;
                return ov;
            }

            offset -= 1;
        }
    }

    OverlapResult ov;
    ov.overlapped = false;
    ov.offset = ov.overlap_len = ov.diff = 0;
    ov.hasGap = false;
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
    string mergedSeq = r1->mSeq->substr(0, len1);
    if(ov.offset > 0) {
        mergedSeq += rr2->mSeq->substr(ol, len2);
    }

    string mergedQual = r1->mQuality->substr(0, len1);
    if(ov.offset > 0) {
        mergedQual += rr2->mQuality->substr(ol, len2);
    }

    delete rr2;

    string name = *(r1->mName) + " merged_" + to_string(len1) + "_" + to_string(len2);
    string strand = *(r1->mStrand);
    if (strand != "+") {
      strand = strand + " merged_" + to_string(len1) + "_" + to_string(len2);
    }
    Read* mergedRead = new Read(new string(name), new string(mergedSeq), new string(strand), new string(mergedQual));

    return mergedRead;
}

bool OverlapAnalysis::test(){
    //Sequence r1("CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGCCGCTGGAGGTCTCCC");
    //Sequence r2("ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGCCCGTAGGCGCGGCTCCC");

    string* r1 = new string("CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGC");
    string* r2 = new string("ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGTCC");
    string* qual1 = new string("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
    string* qual2 = new string("#########################################################################################");
    
    OverlapResult ov = OverlapAnalysis::analyze(r1, r2, 2, 30, 0.2);

    Read read1(new string("name1"), r1, new string("+"), qual1);
    Read read2(new string("name2"), r2, new string("+"), qual2);

    Read* mergedRead = OverlapAnalysis::merge(&read1, &read2, ov);
    mergedRead->print();

    return ov.overlapped && ov.offset == 10 && ov.overlap_len == 79 && ov.diff == 1;
}
