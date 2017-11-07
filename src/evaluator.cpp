#include "evaluator.h"
#include "fastqreader.h"

Evaluator::Evaluator(Options* opt){
    mOptions = opt;
}


Evaluator::~Evaluator(){
}

void Evaluator::evaluate(long& readNum) {
    FastqReader reader(mOptions->in1);

    const long READ_LIMIT = 100000;
    long records = 0;
    long bytes = 0;

    bool reachedEOF = false;
    while(records < READ_LIMIT) {
        Read* r = reader.read();
        if(!r) {
            reachedEOF = true;
            break;
        }
        records++;
    }

    if(reachedEOF){
        readNum = records;
        return ;
    }

    size_t bytesRead;
    size_t bytesTotal;

    reader.getBytes(bytesRead, bytesTotal);

    double bytesPerRead = (double)bytesRead / (double) records;
    readNum = (long) (bytesTotal / bytesPerRead);
}