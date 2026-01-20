#ifndef MATCHER_H
#define MATCHER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

class Matcher{
public:
    Matcher();
    ~Matcher();

    static bool matchWithOneInsertion(const char* insData, const char* normalData, int cmplen, int diffLimit);
    static int diffWithOneInsertion(const char* insData, const char* normalData, int cmplen, int diffLimit);


};


#endif