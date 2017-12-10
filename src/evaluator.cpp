#include "evaluator.h"
#include "fastqreader.h"
#include <map>
#include <memory.h>

Evaluator::Evaluator(Options* opt){
    mOptions = opt;
}


Evaluator::~Evaluator(){
}

bool Evaluator::isTwoColorSystem() {
    FastqReader reader(mOptions->in1);

    Read* r = reader.read();

    if(!r)
        return false;

    // NEXTSEQ500, NEXTSEQ 550, NOVASEQ
    if(starts_with(r->mName, "@NS") || starts_with(r->mName, "@NB") || starts_with(r->mName, "@A0")) {
        delete r;
        return true;
    }

    delete r;
    return false;
}

void Evaluator::evaluateSeqLen() {
    if(!mOptions->in1.empty())
        mOptions->seqLen1 = computeSeqLen(mOptions->in1);
    if(!mOptions->in2.empty())
        mOptions->seqLen2 = computeSeqLen(mOptions->in2);
}

int Evaluator::computeSeqLen(string filename) {
    FastqReader reader(filename);

    long records = 0;
    bool reachedEOF = false;

    // get seqlen
    int seqlen=0;
    while(records < 1000) {
        Read* r = reader.read();
        if(!r) {
            reachedEOF = true;
            break;
        }
        int rlen = r->length();
        if(rlen > seqlen)
            seqlen = rlen;
        records ++;
        delete r;
    }

    return seqlen;
}

void Evaluator::computeOverRepSeq(string filename, map<string, long>& hotseqs, int seqlen) {
    FastqReader reader(filename);

    map<string, long> seqCounts;

    const long BASE_LIMIT = 151 * 10000;
    long records = 0;
    long bases = 0;
    bool reachedEOF = false;

    while(bases < BASE_LIMIT) {
        Read* r = reader.read();
        if(!r) {
            reachedEOF = true;
            break;
        }
        int rlen = r->length();
        bases += rlen;
        records ++;
        // 10, 20, 40, 80, 150

        int steps[5] = {10, 20, 40, 100, seqlen-1};
        
        for(int s=0; s<5; s++) {
            int step = steps[s];
            for(int i=0; i<rlen-step; i++) {
                string seq = r->mSeq.mStr.substr(i, step);
                if(seqCounts.count(seq)>0)
                    seqCounts[seq]++;
                else
                    seqCounts[seq]=1;
            }
        }

        delete r;
    }
    
    map<string, long>::iterator iter;
    for(iter = seqCounts.begin(); iter!=seqCounts.end(); iter++) {
        string seq = iter->first;
        long count = iter->second;

        if(seq.length() >= seqlen-1) {
            if(count >= 3) {
                hotseqs[seq]=count;
            }
        } else if(seq.length() >= 100) {
            if(count >= 5) {
                hotseqs[seq]=count;
            }
        } else if(seq.length() >= 40) {
            if(count >= 20) {
                hotseqs[seq]=count;
            }
        } else if(seq.length() >= 20) {
            if(count >= 100) {
                hotseqs[seq]=count;
            }
        } else if(seq.length() >= 10) {
            if(count >= 500) {
                hotseqs[seq]=count;
            }
        }
    }

    // remove substrings
    map<string, long>::iterator iter2;
    iter = hotseqs.begin(); 
    while(iter!=hotseqs.end()) {
        string seq = iter->first;
        long count = iter->second;
        bool isSubString = false;
        for(iter2 = hotseqs.begin(); iter2!=hotseqs.end(); iter2++) {
            string seq2 = iter2->first;
            long count2 = iter2->second;
            if(seq != seq2 && seq2.find(seq) != string::npos && count / count2 < 10) {
                isSubString = true;
                break;
            }
        }
        if(isSubString) {
            hotseqs.erase(iter++);
        } else {
            iter++;
        }
    }

    // output for test
    /*for(iter = hotseqs.begin(); iter!=hotseqs.end(); iter++) {
        cout << iter->first << ": " << iter->second << endl;
    }*/
}

void Evaluator::evaluateOverRepSeqs() {
    if(!mOptions->in1.empty())
        computeOverRepSeq(mOptions->in1, mOptions->overRepSeqs1, mOptions->seqLen1);
    if(!mOptions->in2.empty())
        computeOverRepSeq(mOptions->in2, mOptions->overRepSeqs2, mOptions->seqLen2);
}

void Evaluator::evaluateReadNum(long& readNum) {
    FastqReader reader(mOptions->in1);

    const long READ_LIMIT = 1024*1024;
    const long BASE_LIMIT = 151 * 1024*1024;
    long records = 0;
    long bases = 0;
    size_t firstReadPos = 0;

    size_t bytesRead;
    size_t bytesTotal;

    bool reachedEOF = false;
    bool first = true;
    while(records < READ_LIMIT && bases < BASE_LIMIT) {
        Read* r = reader.read();
        if(!r) {
            reachedEOF = true;
            break;
        }
        if(first) {
            reader.getBytes(bytesRead, bytesTotal);
            firstReadPos = bytesRead;
            first = false;
        }
        records++;
        bases += r->length();
        delete r;
    }

    if(reachedEOF){
        readNum = records;
        return ;
    }

    reader.getBytes(bytesRead, bytesTotal);

    double bytesPerRead = (double)(bytesRead - firstReadPos) / (double) records;
    // increase it by 1% since the evaluation is usually a bit lower due to bad quality causes lower compression rate
    readNum = (long) (bytesTotal*1.01 / bytesPerRead);
}

string Evaluator::evaluateRead1AdapterAndReadNum(long& readNum) {
    FastqReader reader(mOptions->in1);
    // stat up to 1M reads
    const long READ_LIMIT = 1024*1024;
    const long BASE_LIMIT = 151 * 1024*1024;
    long records = 0;
    long bases = 0;
    size_t firstReadPos = 0;

    size_t bytesRead;
    size_t bytesTotal;

    // we have to shift last cycle for evaluation since it is so noisy, especially for Illumina data
    const int shiftTail = max(1, mOptions->trim.tail1);

    // we count the last [2, 9] bp of each read
    // why we add trim_tail here? since the last cycle are usually with low quality and should be trimmed
    const int keylen = 10;
    int size = 1 << (keylen*2 );
    unsigned int* counts = new unsigned int[size];
    memset(counts, 0, sizeof(unsigned int)*size);
    bool reachedEOF = false;
    bool first = true;
    while(records < READ_LIMIT && bases < BASE_LIMIT) {
        Read* r = reader.read();
        if(!r) {
            reachedEOF = true;
            break;
        }
        if(first) {
            reader.getBytes(bytesRead, bytesTotal);
            firstReadPos = bytesRead;
            first = false;
        }
        int rlen = r->length();
        bases += rlen;
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

    if(reachedEOF){
        readNum = records;
        return finalAdapter;
    }

    // by the way, update readNum so we don't need to evaluate it if splitting output is enabled
    reader.getBytes(bytesRead, bytesTotal);

    double bytesPerRead = (double)(bytesRead - firstReadPos) / (double) records;
    // increase it by 1% since the evaluation is usually a bit lower due to bad quality causes lower compression rate
    readNum = (long) (bytesTotal*1.01 / bytesPerRead);

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