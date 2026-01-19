#include "evaluator.h"
#include "fastqreader.h"
#include <map>
#include <memory.h>
#include "nucleotidetree.h"
#include "knownadapters.h"

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

    // Two-color system instrument prefixes
    // See https://github.com/OpenGene/fastp/pull/508 for a detailed discussion
    const vector<string> twoColorPrefixes = {
        "@FS",  // iSeq 100
        "@MN", "@SH",  // MiniSeq
        "@NS", "@NB",  // NextSeq 500/550
        "@NDX",  // NextSeq 550DX
        "@VL", "@VH",  // NextSeq 1000/2000
        "@A", "@NA",  // NovaSeq 6000
        "@LH"  // NovaSeq X
    };

    for (const string& prefix : twoColorPrefixes) {
        if (starts_with(r->mName, prefix)) {
            delete r;
            return true;
        }
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

        int steps[5] = {10, 20, 40, 100, min(150,seqlen-2)};
        
        for(int s=0; s<5; s++) {
            int step = steps[s];
            for(int i=0; i<rlen-step; i++) {
                string seq = r->mSeq->substr(i, step);
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
        cerr << iter->first << ": " << iter->second << endl;
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

    const long READ_LIMIT = 512*1024;
    const long BASE_LIMIT = 151 * 512*1024;
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

    readNum = 0;
    if(reachedEOF){
        readNum = records;
    } else if(records>0) {
        // by the way, update readNum so we don't need to evaluate it if splitting output is enabled
        reader.getBytes(bytesRead, bytesTotal);
        double bytesPerRead = (double)(bytesRead - firstReadPos) / (double) records;
        // increase it by 1% since the evaluation is usually a bit lower due to bad quality causes lower compression rate
        readNum = (long) (bytesTotal*1.01 / bytesPerRead);
    }
}

string Evaluator::checkKnownAdapters(Read** reads, long num) {
    map<string, string> knownAdapters = getKnownAdapter();
    map<string, int> possibleCounts;
    map<string, int> mismatches;

    // for performance, up to 100k reads and 100M bases
    const int MAX_CHECK_READS = 100000;
    const int MAX_CHECK_BASES = MAX_CHECK_READS*1000;
    // if we hit for 1000 times, then we exit
    const int MAX_HIT = 1000;

    const int matchReq = 8;
    const int allowOneMismatchForEach = 16;

    map<string, string>::iterator iter;

    for(iter = knownAdapters.begin(); iter!= knownAdapters.end(); iter++) {
        string adapter = iter->first;
        possibleCounts[adapter] = 0;
        mismatches[adapter] = 0;
    }

    long checkedReads = 0;
    long checkedBases = 0;
    int curMaxCount = 0;
    for (long i=0; i<num; i++) {
        Read* r = reads[i];
        const char* rdata = r->mSeq->c_str();
        int rlen = r->length();

        checkedReads++;
        checkedBases+= rlen;
        if(checkedReads > MAX_CHECK_READS || checkedBases > MAX_CHECK_BASES)
            break;
        if(curMaxCount > MAX_HIT)
            break;
        for(iter = knownAdapters.begin(); iter!= knownAdapters.end(); iter++) {
            string adapter = iter->first;
            const char* adata = adapter.c_str();
            int alen = adapter.length();
            if(alen >= rlen)
                continue;
            // this one is not the candidate, skip it for speedup
            if(curMaxCount > 20 && possibleCounts[adapter] <curMaxCount/10) {
                continue; 
            }
            for(int pos = 0; pos<rlen-matchReq; pos++) {
                int cmplen = min(rlen - pos, alen);
                int allowedMismatch = cmplen/allowOneMismatchForEach;
                int mismatch = 0;
                bool matched = true;
                for(int i=0; i<cmplen; i++) {
                    if( adata[i] != rdata[i+pos] ){
                        mismatch++;
                        if(mismatch > allowedMismatch) {
                            matched = false;
                            break;
                        }
                    }
                }
                if(matched) {
                    possibleCounts[adapter]++;
                    if(curMaxCount < possibleCounts[adapter])
                        curMaxCount = possibleCounts[adapter];
                    mismatches[adapter] += mismatch;
                    break;
                }
            }
        }
    }

    string adapter = "";
    int maxCount = 0;
    map<string, int>::iterator iter2;
    for(iter2 = possibleCounts.begin(); iter2 != possibleCounts.end(); iter2++) {
        if(iter2->second > maxCount) {
            adapter = iter2->first;
            maxCount = iter2->second;
        }
    }
    if(maxCount > checkedReads/50 || (maxCount > checkedReads/200 && mismatches[adapter] < checkedReads)) {
        cerr << knownAdapters[adapter] << endl;
        cerr << adapter << endl;
        return adapter;
    }
    return "";
}

string Evaluator::evalAdapterAndReadNum(long& readNum, bool isR2) {
    string filename = mOptions->in1;
    if(isR2)
        filename = mOptions->in2;
    FastqReader reader(filename);
    // stat up to 256K reads
    const long READ_LIMIT = 256*1024;
    const long BASE_LIMIT = 151 * READ_LIMIT;
    long records = 0;
    long bases = 0;
    size_t firstReadPos = 0;

    size_t bytesRead;
    size_t bytesTotal;

    Read** loadedReads = new Read*[READ_LIMIT];
    memset(loadedReads, 0, sizeof(Read*)*READ_LIMIT);
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
        loadedReads[records] = r;
        records++;
    }

    readNum = 0;
    if(reachedEOF){
        readNum = records;
    } else if(records>0) {
        // by the way, update readNum so we don't need to evaluate it if splitting output is enabled
        reader.getBytes(bytesRead, bytesTotal);
        double bytesPerRead = (double)(bytesRead - firstReadPos) / (double) records;
        // increase it by 1% since the evaluation is usually a bit lower due to bad quality causes lower compression rate
        readNum = (long) (bytesTotal*1.01 / bytesPerRead);
    }

    // we need at least 10000 valid records to evaluate
    if(records < 10000) {
        for(int r=0; r<records; r++) {
            delete loadedReads[r];
            loadedReads[r] = NULL;
        }
        delete[] loadedReads;
        return "";
    }

    string knownAdapter = checkKnownAdapters(loadedReads, records);
     if(knownAdapter.size() > 8) {
        for(int r=0; r<records; r++) {
            delete loadedReads[r];
            loadedReads[r] = NULL;
        }
        delete[] loadedReads;
        return knownAdapter;
    }

    // we have to shift last cycle for evaluation since it is so noisy, especially for Illumina data
    const int shiftTail = max(1, mOptions->trim.tail1);

    // why we add trim_tail here? since the last cycle are usually with low quality and should be trimmed
    const int keylen = 10;
    int size = 1 << (keylen*2 );
    unsigned int* counts = new unsigned int[size];
    memset(counts, 0, sizeof(unsigned int)*size);
    for(int i=0; i<records; i++) {
        Read* r = loadedReads[i];
        const char* data = r->mSeq->c_str();
        int key = -1;
        for(int pos = 20; pos <= r->length()-keylen-shiftTail; pos++) {
            key = seq2int(r->mSeq, pos, keylen, key);
            if(key >= 0) {
                counts[key]++;
            }
        }
    }

    // set AAAAAAAAAA = 0;
    counts[0] = 0;

    // get the top N
    const int topnum = 10;
    int topkeys[topnum] = {0};
    long total = 0;
    for(int k=0; k<size; k++) {
        int atcg[4] = {0};
        for(int i=0; i<keylen; i++) {
            int baseOfBit = (k >> (i*2)) & 0x03;
            atcg[baseOfBit]++;
        }
        bool lowComplexity = false;
        for(int b=0; b<4; b++) {
            if(atcg[b] >= keylen-4)
                lowComplexity=true;
        }
        if(lowComplexity)
            continue;
        // too many GC
        if(atcg[2] + atcg[3] >= keylen-2)
            continue;

        // starts with GGGG
        if( k>>12 == 0xff)
            continue;

        unsigned int val = counts[k];
        total += val;
        for(int t=topnum-1; t>=0; t--) {
            // reach the middle
            if(val < counts[topkeys[t]]){
                if(t<topnum-1) {
                    for(int m=topnum-1; m>t+1; m--) {
                        topkeys[m] = topkeys[m-1];
                    }
                    topkeys[t+1] = k;
                }
                break;
            } else if(t == 0) { // reach the top
                for(int m=topnum-1; m>t; m--) {
                    topkeys[m] = topkeys[m-1];
                }
                topkeys[t] = k;
            }
        }
    }

    const int FOLD_THRESHOLD = 20;
    for(int t=0; t<topnum; t++) {
        int key = topkeys[t];
        string seq = int2seq(key, keylen);
        if(key == 0)
            continue;
        long count = counts[key];
        if(count<10 || count*size < total * FOLD_THRESHOLD)
            break;
        // skip low complexity seq
        int diff = 0;
        for(int s=0; s<seq.length() - 1; s++) {
            if(seq[s] != seq[s+1])
                diff++;
        }
        if(diff <3){
            continue;
        }
        string adapter = getAdapterWithSeed(key, loadedReads, records, keylen);
        if(!adapter.empty()){
            delete[] counts;
            for(int r=0; r<records; r++) {
                delete loadedReads[r];
                loadedReads[r] = NULL;
            }
            delete[] loadedReads;
            return adapter;
        }
    }

    delete[] counts;
    for(int r=0; r<records; r++) {
        delete loadedReads[r];
        loadedReads[r] = NULL;
    }
    delete[] loadedReads;
    return "";

}

string Evaluator::getAdapterWithSeed(int seed, Read** loadedReads, long records, int keylen) {
    // we have to shift last cycle for evaluation since it is so noisy, especially for Illumina data
    const int shiftTail = max(1, mOptions->trim.tail1);
    const int MAX_SEARCH_LENGTH = 500;
    NucleotideTree forwardTree(mOptions);
    // forward search
    for(int i=0; i<records; i++) {
        Read* r = loadedReads[i];
        const char* data = r->mSeq->c_str();
        int key = -1;
        for(int pos = 20; pos <= r->length()-keylen-shiftTail && pos <MAX_SEARCH_LENGTH; pos++) {
            key = seq2int(r->mSeq, pos, keylen, key);
            if(key == seed) {
                forwardTree.addSeq(r->mSeq->substr(pos+keylen, r->length()-keylen-shiftTail-pos));
            }
        }
    }
    bool reachedLeaf = true;
    string forwardPath = forwardTree.getDominantPath(reachedLeaf);

    NucleotideTree backwardTree(mOptions);
    // backward search
    for(int i=0; i<records; i++) {
        Read* r = loadedReads[i];
        const char* data = r->mSeq->c_str();
        int key = -1;
        for(int pos = 20; pos <= r->length()-keylen-shiftTail && pos <MAX_SEARCH_LENGTH; pos++) {
            key = seq2int(r->mSeq, pos, keylen, key);
            if(key == seed) {
                string seq =  r->mSeq->substr(0, pos);
                string rcseq = reverse(seq);
                backwardTree.addSeq(rcseq);
            }
        }
    }
    string backwardPath = backwardTree.getDominantPath(reachedLeaf);

    string adapter = reverse(backwardPath) + int2seq(seed, keylen) + forwardPath;
    if(adapter.length()>60)
        adapter.resize(60);

    string matchedAdapter = matchKnownAdapter(adapter);
    if(!matchedAdapter.empty()) {
        map<string, string> knownAdapters = getKnownAdapter();
        cerr << knownAdapters[matchedAdapter] << endl << matchedAdapter << endl;
        return matchedAdapter;
    } else {
        if(reachedLeaf) {
            cerr << adapter << endl;
            return adapter;
        } else {
            return "";
        }
    }
}

string Evaluator::matchKnownAdapter(string seq) {
    map<string, string> knownAdapters = getKnownAdapter();
    map<string, string>::iterator iter;
    for(iter = knownAdapters.begin(); iter != knownAdapters.end(); iter++) {
        string adapter = iter->first;
        string desc = iter->second;
        if(seq.length()<adapter.length()) {
            continue;
        }
        int diff = 0;
        for(int i=0; i<adapter.length() && i<seq.length(); i++) {
            if(adapter[i] != seq[i])
                diff++;
        }
        if(diff == 0)
            return adapter;
    }
    return "";
}

string Evaluator::int2seq(unsigned int val, int seqlen) {
    char bases[4] = {'A', 'T', 'C', 'G'};
    string ret(seqlen, 'N');
    int done = 0;
    while(done < seqlen) {
        ret[seqlen - done - 1] = bases[val & 0x03];
        val = (val >> 2);
        done++;
    }
    return ret;
}

int Evaluator::seq2int(string* seq, int pos, int keylen, int lastVal) {
    return seq2int(*seq, pos, keylen, lastVal);
}

int Evaluator::seq2int(string& seq, int pos, int keylen, int lastVal) {
    int rlen = seq.length();
    if(lastVal >= 0) {
        const int mask = (1 << (keylen*2 )) - 1;
        int key = (lastVal<<2) & mask;
        char base = seq[pos + keylen - 1];
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
                return -1;
        }
        return key;
    } else {
        int key = 0;
        for(int i=pos; i<keylen+pos; i++) {
            key = (key << 2);
            char base = seq[i];
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
                    return -1;
            }
        }
        return key;
    }
}

bool Evaluator::test() {
    Evaluator eval(NULL);
    string s = "ATCGATCGAT";
    cerr << eval.int2seq(eval.seq2int(s, 0, 10, -1), 10) << endl;
    return eval.int2seq(eval.seq2int(s, 0, 10, -1), 10) == s;
}
