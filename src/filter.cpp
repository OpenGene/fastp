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
    if(r == NULL) {
        return FAIL_LENGTH;
    }

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

Read* Filter::trimAndCut(Read* r, int front, int tail) {
    // return the same read for speed if no change needed
    if(front == 0 && tail == 0 && !mOptions->qualityCut.enabled5 && !mOptions->qualityCut.enabled3)
        return r;


    int rlen = r->length() - front - tail ; 
    if (rlen < 0)
        return NULL;

    if(front == 0 && !mOptions->qualityCut.enabled5 && !mOptions->qualityCut.enabled3){
        r->resize(rlen);
        return r;
    } else if(!mOptions->qualityCut.enabled5 && !mOptions->qualityCut.enabled3){
        r->mSeq.mStr = r->mSeq.mStr.substr(front, rlen);
        r->mQuality = r->mQuality.substr(front, rlen);
        return r;
    }

    // need quality cutting

    int w = mOptions->qualityCut.windowSize;
    int l = r->length();
    const char* qualstr = r->mQuality.c_str();
    const char* seq = r->mSeq.mStr.c_str();
    // quality cutting forward
    if(mOptions->qualityCut.enabled5) {
        int s = front;
        if(l - front - tail - w <= 0)
            return NULL;

        int totalQual = 0;

        // preparing rolling
        for(int i=0; i<w-1; i++)
            totalQual += qualstr[s+i];

        for(s=front; s+w<l-tail; s++) {
            totalQual += qualstr[s+w-1];
            // rolling
            if(s > front) {
                totalQual -= qualstr[s-1];
            }
            // add 33 for phred33 transforming
            if((double)totalQual / (double)w >= 33 + mOptions->qualityCut.quality)
                break;
        }

        // the trimming in front is forwarded and rlen is recalculated
        if(s >0 )
            s = s+w-1;
        while(s<l && seq[s] == 'N')
            s++;
        front = s;
        rlen = l - front - tail;
    }

    // quality cutting backward
    if(mOptions->qualityCut.enabled3) {
        if(l - front - tail - w <= 0)
            return NULL;

        int totalQual = 0;
        int t = l - tail - 1;

        // preparing rolling
        for(int i=0; i<w-1; i++)
            totalQual += qualstr[t-i];

        for(t=l-tail-1; t-w>=front; t--) {
            totalQual += qualstr[t-w+1];
            // rolling
            if(t < l-tail-1) {
                totalQual -= qualstr[t+1];
            }
            // add 33 for phred33 transforming
            if((double)totalQual / (double)w >= 33 + mOptions->qualityCut.quality)
                break;
        }

        if(t < l-1)
            t = t-w+1;
        while(t>=0 && seq[t] == 'N')
            t--;
        rlen = t - front + 1;
    }

    if(rlen <= 0 || front >= l-1)
        return NULL;

    r->mSeq.mStr = r->mSeq.mStr.substr(front, rlen);
    r->mQuality = r->mQuality.substr(front, rlen);

    return r;
}

bool Filter::test() {
    Read r("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTT",
        "+",
        "/////CCCCCCCCCCCC////CCCCCCCCCCCCCC////E");
    Options opt;
    opt.qualityCut.enabled5 = true;
    opt.qualityCut.enabled3 = true;
    opt.qualityCut.windowSize = 4;
    opt.qualityCut.quality = 20;
    Filter filter(&opt);
    Read* ret = filter.trimAndCut(&r, 0, 1);
    
    return ret->mSeq.mStr == "TAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAAT"
        && ret->mQuality == "//CCCCCCCCCCCC////CCCCCCCCCCCCCC//";
}