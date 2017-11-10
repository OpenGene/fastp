#ifndef EVALUATOR_H
#define EVALUATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "options.h"

using namespace std;

class Evaluator{
public:
    Evaluator(Options* opt);
    ~Evaluator();
    // evaluate how many reads are stored in the input file
    void evaluateReads(long& readNum);
    string evaluateRead1Adapter();
private:
    Options* mOptions;
    string int2seq(unsigned int val, int seqlen);
    unsigned int seq2int(string& seq, int seqlen, bool& valid);
};


#endif