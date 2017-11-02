#include "processor.h"
#include "peprocessor.h"
#include "seprocessor.h"

Filter::Filter(Options* opt){
    mOptions = opt;
}


Filter::~Filter(){
}

bool Filter::passFilter(Read* r, int lowQualNum, int nBaseNum) {
    if(mOptions->qualfilter.enabled) {
        if(lowQualNum > mOptions->qualfilter.unqualifiedBaseLimit )
            return false;
        else if(nBaseNum > mOptions->qualfilter.nBaseLimit )
            return false;
    }

    return true;
}