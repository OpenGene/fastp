#ifndef FILTER_H
#define FILTER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "options.h"

using namespace std;

class Filter{
public:
    Filter(Options* opt);
    ~Filter();
    bool passFilter(Read* r, int lowQualNum, int nBaseNum);

private:
    Options* mOptions;
};


#endif