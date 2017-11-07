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
    void evaluate(long& readNum);
private:
    Options* mOptions;
};


#endif