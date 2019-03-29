#ifndef FILTER_H
#define FILTER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "options.h"
#include "read.h"

using namespace std;

class Filter{
public:
    Filter(Options* opt);
    ~Filter();
    int passFilter(Read* r);
    bool passLowComplexityFilter(Read* r);
    Read* trimAndCut(Read* r, int front, int tail, int& frontTrimmed);
    bool filterByIndex(Read* r);
    bool filterByIndex(Read* r1, Read* r2);
    static bool test();

private:
    bool match(vector<string>& list, string target, int threshold);

private:
    Options* mOptions;
};


#endif