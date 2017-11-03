#include "processor.h"
#include "peprocessor.h"
#include "seprocessor.h"
#include "overlapanalysis.h"

Filter::Filter(Options* opt){
    mOptions = opt;
}


Filter::~Filter(){
}

int Filter::passFilter(Read* r, int lowQualNum, int nBaseNum) {

    if(mOptions->qualfilter.enabled) {
        if(lowQualNum > (mOptions->qualfilter.unqualifiedPercentLimit * r->length() / 100.0) )
            return FAIL_QUALITY;
        else if(nBaseNum > mOptions->qualfilter.nBaseLimit )
            return FAIL_N_BASE;
    }

    if(mOptions->lengthFilter.enabled) {
        if(r->length() < mOptions->lengthFilter.requiredLength)
            return FAIL_LENGTH;
    }

    return PASS_FILTER;
}

Read* Filter::trimAndCutAdapter(Read* r) {
    // return the same read for speed if no trimming applied
    if(mOptions->trim.front == 0 && mOptions->trim.tail == 0)
        return r;

    int rlen = r->length() - mOptions->trim.front - mOptions->trim.tail ; 
    if (rlen < 0)
        return NULL;

    Read* ret = new Read(r->mName, r->mSeq.mStr.substr(mOptions->trim.front, rlen), r->mStrand, r->mQuality.substr(mOptions->trim.front, rlen));
    return ret;
}