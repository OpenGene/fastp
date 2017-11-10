#include "evaluator.h"
#include "fastqreader.h"
#include <map>

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
    // stat up to 1M reads
    const long READ_LIMIT = 1024*1024;
    long records = 0;
    // we have to shift last cycle for evaluation since it is so noisy, especially for Illumina data
    const int shiftTail = max(1, mOptions->trim.tail1);

    // we count the last [2, 9] bp of each read
    // why we add trim_tail here? since the last cycle are usually with low quality and should be trimmed
    const int keylen = 10;
    int size = 1 << (keylen*2 );
    unsigned int* counts = new unsigned int[size];
    memset(counts, 0, sizeof(unsigned int)*size);
    bool reachedEOF = false;
    while(records < READ_LIMIT) {
        Read* r = reader.read();
        if(!r) {
            reachedEOF = true;
            break;
        }
        int rlen = r->length();
        if(rlen < keylen + 1 + shiftTail)
            continue;

        const char* data = r->mSeq.mStr.c_str();
        bool valid = true;
        unsigned int key = 0;
        for(int i=0; i<keylen; i++) {
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
        delete r;
    }

    // we need at least 10000 valid records to evaluate
    if(records < 10000) {
        delete counts;
        return "";
    }

    int repeatReq = 0.0001 * records;

    // initialize candidates
    map<string, unsigned int> candidates;
    for(int i=0; i<size; i++) {
        if(counts[i] >= repeatReq) {
            string seq = int2seq(i, keylen);
            // remove low complexity seq
            int diff = 0;
            for(int s=0; s<seq.length() - 1; s++) {
                if(seq[s] != seq[s+1])
                    diff++;
            }
            if(diff >=2){
                candidates[seq] = counts[i];
                //cout << seq << ": " << candidates[seq] << endl;
            }
        }
    }

    map<string, unsigned int>::iterator iter;

    // remove the fake ones have only first base different
    vector<string> needToDelete;
    for(iter = candidates.begin(); iter!=candidates.end(); iter++) {
        string seq = iter->first;
        char bases[4] = {'A', 'T', 'C', 'G'};
        int num = 0;
        for(int b=0; b<4; b++) {
            seq[0] = bases[b];
            if(candidates.count(seq) > 0)
                num++;
        }
        if(num >=2 ) {
            needToDelete.push_back(iter->first);
        }
    }
    for(int i=0; i<needToDelete.size(); i++) {
        candidates.erase(needToDelete[i]);
    }

    map<string, unsigned int>::iterator iter1;
    map<string, unsigned int>::iterator iter2;

    while(true) {
        bool changed = false;
        for(iter1 = candidates.begin(); iter1!=candidates.end(); iter1++) {
            bool aligned = false;
            for(iter2 = candidates.begin(); iter2!=candidates.end(); iter2++) {
                if(iter1 == iter2)
                    continue;

                string a1 = iter1->first;
                string a2 = iter2->first;
                int len1 = a1.length();
                int len2 = a2.length();
                int overlap = keylen - 1;
                //cout << a1 << ":" << a2 << endl;

                // check identidal
                bool identical = true;
                for(int o=0; o<overlap; o++) {
                    identical &= (a1[len1 - overlap + o] == a2[o]);
                }

                if(identical) {
                    // merge them
                    string mergedAdapter = a1 + a2.substr(overlap, len2-overlap);
                    int mergedCount = iter1->second + iter2->second;
                    candidates.erase(a1);
                    candidates.erase(a2);
                    candidates[mergedAdapter] = mergedCount;
                    aligned = true;
                    break;
                }

            }
            if(aligned) {
                changed = true;
                break;
            }
        }
        if(changed == false)
            break;
    }

    // find the longest adapter
    int largest = 0;
    string finalAdapter = "";
    for(iter = candidates.begin(); iter!=candidates.end(); iter++) {
        if(iter->second > largest) {
            largest = iter->second;
            finalAdapter = iter->first;
        }
    }

    delete counts;
    return finalAdapter;

}

string Evaluator::int2seq(unsigned int val, int seqlen) {
    char bases[4] = {'A', 'T', 'C', 'G'};
    string ret(seqlen, 'N');
    int done = 0;
    while(done < seqlen) {
        ret[done] = bases[val & 0x03];
        val = (val >> 2);
        done++;
    }
    return ret;
}

unsigned int Evaluator::seq2int(string& seq, int keylen, bool& valid) {
    valid = true;
    unsigned int key = 0;
    int rlen = seq.length();
    for(int i=0; i<keylen; i++) {
        key = (key << 2);
        char base = seq[rlen - i - 1];
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
                return 0;
        }
    }
    return key;
}