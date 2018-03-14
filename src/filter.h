#ifndef FILTER_H
#define FILTER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "options.h"
#include "read.h"

using namespace std;

class Filter{
public:
    Filter(Options* opt);
    ~Filter();
    int passFilter(Read* r, int lowQualNum, int nBaseNum);
    bool passLowComplexityFilter(Read* r);
    Read* trimAndCut(Read* r, int front, int tail);
    static bool test();

private:
    Options* mOptions;
};


#endif