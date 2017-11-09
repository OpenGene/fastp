#include "evaluator.h"
#include "fastqreader.h"

Evaluator::Evaluator(Options* opt){
    mOptions = opt;
}


Evaluator::~Evaluator(){
}

void Evaluator::evaluateReads(long& readNum) {
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
        delete r;
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

string Evaluator::evaluateRead1Adapter() {
    FastqReader reader(mOptions->in1);
    // stat up to 64K reads
    const long READ_LIMIT = 65536;
    long records = 0;
    const int shiftTail = 1;

    // we count the last [2, 9] bp of each read
    // why we add trim_tail here? since the last cycle are usually with low quality and should be trimmed
    unsigned short* counts = new unsigned short[65536];
    bool reachedEOF = false;
    while(records < READ_LIMIT) {
        Read* r = reader.read();
        if(!r) {
            reachedEOF = true;
            break;
        }
        int rlen = r->length();
        if(rlen < 9)
            continue;

        const char* data = r->mSeq.mStr.c_str();
        bool valid = true;
        unsigned short key = 0;
        for(int i=0; i<8; i++) {
            key = (key << 2);
            char base = data[rlen - i - 1 - shiftTail];
            switch (base) {
                case 'A':
                    key += 0;
                    break;
                case 'T':
                    key += 1;
                    break;
                case 'C':
                    key += 2;
                    break;
                case 'G':
                    key += 3;
                    break;
                default:
                    // N or anything else
                    valid = false;
                    break;
            }
            if(!valid)
                break;
        }
        if(valid) {
            counts[key]++;
            records++;
        }
    }

    // we need at least 10000 valid records to evaluate
    if(records < 10000)
        return "";

    return "";

    // get those frequency >= 0.0001;
}