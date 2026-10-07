#include "matcher.h"
#include <cstdlib>
#include <string>
#include <iostream>

Matcher::Matcher(){
}


Matcher::~Matcher(){
}

// insData has one base more than normalData. Returns true if, for some split i in
// [1, cmplen-1] (insData[i] being the inserted base),
//     mismatches(insData[0, i), normalData[0, i)) + mismatches(insData(i, cmplen], normalData[i, cmplen))
// is at most diffLimit. Same result as matchWithOneInsertionReference, which builds both
// prefix/suffix mismatch arrays in full; here each scan stops as soon as it can no longer
// be part of a match. On unrelated sequence both stop after a few bases, so a call is O(diffLimit)
// instead of O(cmplen).
bool Matcher::matchWithOneInsertion(const char* insData, const char* normalData, int cmplen, int diffLimit) {
    if (cmplen < 2 || diffLimit < 0)
        return false;
    if (cmplen > MAX_FAST_CMPLEN)
        return matchWithOneInsertionReference(insData, normalData, cmplen, diffLimit);

    const int head = insData[0] != normalData[0] ? 1 : 0;              // left[0]
    const int tail = insData[cmplen] != normalData[cmplen - 1] ? 1 : 0; // right[cmplen-1]

    // left[k] = mismatches over [0, k]. Split i needs left[i-1] + right[i] <= diffLimit, and
    // right[i] >= tail, so only splits with left[i-1] + tail <= diffLimit can match: i <= hi.
    int left[MAX_FAST_CMPLEN];
    int hi = 0;
    int acc = 0;
    for (int k = 0; k < cmplen - 1; k++) {
        acc += insData[k] != normalData[k] ? 1 : 0;
        if (acc + tail > diffLimit)
            break;
        left[k] = acc;
        hi = k + 1;
    }
    if (hi == 0)
        return false;

    // right[j] = mismatches of insData[j+1..cmplen] vs normalData[j..cmplen-1]. left[i-1] >= head,
    // so only splits with right[i] + head <= diffLimit can match: i >= lo.
    int right[MAX_FAST_CMPLEN];
    int lo = cmplen;
    acc = 0;
    for (int j = cmplen - 1; j >= 1; j--) {
        acc += insData[j + 1] != normalData[j] ? 1 : 0;
        if (acc + head > diffLimit)
            break;
        right[j] = acc;
        lo = j;
    }

    for (int i = lo; i <= hi; i++) {
        if (left[i - 1] + right[i] <= diffLimit)
            return true;
    }
    return false;
}

bool Matcher::matchWithOneInsertionReference(const char* insData, const char* normalData, int cmplen, int diffLimit) {
    // accumlated mismatches from left/right
    int accMismatchFromLeft[cmplen];
    int accMismatchFromRight[cmplen];

    // accMismatchFromLeft[0]: head vs. head
    // accMismatchFromRight[cmplen-1]: tail vs. tail
    accMismatchFromLeft[0] = insData[0] == normalData[0] ? 0 : 1;
    accMismatchFromRight[cmplen-1] = insData[cmplen] == normalData[cmplen-1] ? 0 : 1;
    for(int i=1; i<cmplen; i++) {
        if(insData[i] != normalData[i])
            accMismatchFromLeft[i] = accMismatchFromLeft[i-1]+1;
        else
            accMismatchFromLeft[i] = accMismatchFromLeft[i-1];
        
        if(accMismatchFromLeft[i] + accMismatchFromRight[cmplen-1] >diffLimit)
            break;
    }
    for(int i=cmplen - 2; i>=0; i--) {
        if(insData[i+1] != normalData[i])
            accMismatchFromRight[i] = accMismatchFromRight[i+1]+1;
        else
            accMismatchFromRight[i] = accMismatchFromRight[i+1];
        if(accMismatchFromRight[i] + accMismatchFromLeft[0]> diffLimit) {
            for(int p=0; p<i; p++)
                accMismatchFromRight[p] = diffLimit+1;
            break;
        }
    }

    //    insData:     XXXXXXXXXXXXXXXXXXXXXXX[i]XXXXXXXXXXXXXXXXXXXXXXXX
    // normalData:     YYYYYYYYYYYYYYYYYYYYYYY   YYYYYYYYYYYYYYYYYYYYYYYY
    //       diff:    accMismatchFromLeft[i-1] + accMismatchFromRight[i]

    // insertion can be from pos = 1 to cmplen - 1
    for(int i=1; i<cmplen; i++) {
        if(accMismatchFromLeft[i-1] + accMismatchFromRight[cmplen-1]> diffLimit)
            return false;
        int diff = accMismatchFromLeft[i-1] + accMismatchFromRight[i];
        if(diff <= diffLimit)
            return true;
    }

    return false;
}

int Matcher::diffWithOneInsertion(const char* insData, const char* normalData, int cmplen, int diffLimit) {
    // accumlated mismatches from left/right
    int accMismatchFromLeft[cmplen];
    int accMismatchFromRight[cmplen];

    // accMismatchFromLeft[0]: head vs. head
    // accMismatchFromRight[cmplen-1]: tail vs. tail
    accMismatchFromLeft[0] = insData[0] == normalData[0] ? 0 : 1;
    accMismatchFromRight[cmplen-1] = insData[cmplen] == normalData[cmplen-1] ? 0 : 1;
    for(int i=1; i<cmplen; i++) {
        if(insData[i] != normalData[i])
            accMismatchFromLeft[i] = accMismatchFromLeft[i-1]+1;
        else
            accMismatchFromLeft[i] = accMismatchFromLeft[i-1];
        
        if(accMismatchFromLeft[i] + accMismatchFromRight[cmplen-1] >diffLimit)
            break;
    }
    for(int i=cmplen - 2; i>=0; i--) {
        if(insData[i+1] != normalData[i])
            accMismatchFromRight[i] = accMismatchFromRight[i+1]+1;
        else
            accMismatchFromRight[i] = accMismatchFromRight[i+1];
        if(accMismatchFromRight[i] + accMismatchFromLeft[0]> diffLimit) {
            for(int p=0; p<i; p++)
                accMismatchFromRight[p] = diffLimit+1;
            break;
        }
    }

    //    insData:     XXXXXXXXXXXXXXXXXXXXXXX[i]XXXXXXXXXXXXXXXXXXXXXXXX
    // normalData:     YYYYYYYYYYYYYYYYYYYYYYY   YYYYYYYYYYYYYYYYYYYYYYYY
    //       diff:    accMismatchFromLeft[i-1] + accMismatchFromRight[i]

    int minDiff = 100000000;
    // insertion can be from pos = 1 to cmplen - 1
    for(int i=1; i<cmplen; i++) {
        if(accMismatchFromLeft[i-1] + accMismatchFromRight[cmplen-1]> diffLimit)
            return -1; // -1 means higher than diffLimit
        int diff = accMismatchFromLeft[i-1] + accMismatchFromRight[i];
        if(diff <= minDiff)
            minDiff = diff;
    }

    return minDiff;
}
bool Matcher::test() {
    // matchWithOneInsertion must agree with the reference implementation everywhere:
    // random pairs, pairs that are one insertion (plus a few substitutions) apart, and limits
    // from below zero up to cmplen.
    srand(721);
    const char* bases = "ACGTN";
    for (int iter = 0; iter < 200000; iter++) {
        int cmplen = 1 + rand() % 70;
        std::string normal(cmplen, 'A');
        for (int k = 0; k < cmplen; k++)
            normal[k] = bases[rand() % 5];
        std::string ins;
        if (iter % 2 == 0) {
            ins.resize(cmplen + 1);
            for (int k = 0; k <= cmplen; k++)
                ins[k] = bases[rand() % 5];
        } else {
            int at = 1 + rand() % cmplen;
            ins = normal.substr(0, at) + bases[rand() % 4] + normal.substr(at);
            int subs = rand() % 4;
            for (int k = 0; k < subs; k++)
                ins[rand() % ins.length()] = bases[rand() % 5];
        }
        int limit = rand() % (cmplen / 4 + 3) - 1;
        bool fast = matchWithOneInsertion(ins.c_str(), normal.c_str(), cmplen, limit);
        bool ref = matchWithOneInsertionReference(ins.c_str(), normal.c_str(), cmplen, limit);
        if (fast != ref) {
            std::cerr << "matchWithOneInsertion mismatch: ins=" << ins << " normal=" << normal
                      << " cmplen=" << cmplen << " limit=" << limit << " fast=" << fast << " ref=" << ref << std::endl;
            return false;
        }
    }
    return true;
}
