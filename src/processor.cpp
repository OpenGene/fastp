#include "processor.h"
#include "peprocessor.h"
#include "seprocessor.h"

Processor::Processor(Options* opt){
    mOptions = opt;
}


Processor::~Processor(){
}

bool Processor::process() {
    if(mOptions->isPaired()) {
        PairEndProcessor p(mOptions);
        p.process();
    } else {
        SingleEndProcessor p(mOptions);
        p.process();
    }

    return true;
}