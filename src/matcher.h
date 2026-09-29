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
    // straightforward version kept to check the fast one (Matcher::test) and for cmplen > MAX_FAST_CMPLEN
    static bool matchWithOneInsertionReference(const char* insData, const char* normalData, int cmplen, int diffLimit);
    static int diffWithOneInsertion(const char* insData, const char* normalData, int cmplen, int diffLimit);
    static bool test();

    static const int MAX_FAST_CMPLEN = 512;


};


#endif