#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>

using namespace std;

class Sequence{
public:
    Sequence();
    Sequence(string*  seq);
    ~Sequence();
    void print();
    int length();
    Sequence reverseComplement();

    Sequence operator~();

    static bool test();
    static string reverseComplement(string* origin);

public:
    string*  mStr;
};

#endif