#include "matcher.h"

Matcher::Matcher(){
}


Matcher::~Matcher(){
}

bool Matcher::matchWithOneInsertion(const char* insData, const char* normalData, int cmplen, int diffLimit) {
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