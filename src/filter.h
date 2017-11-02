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
    bool passFilter(Read* r, int lowQualNum, int nBaseNum);
    Read* trimAndCutAdapter(Read* r);

private:
    Options* mOptions;
};


#endif