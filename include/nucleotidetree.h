#ifndef NUCLEICTREE_H
#define NUCLEICTREE_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <memory.h>
#include "options.h"

using namespace std;

// (A,T,C,G,N) & 0X07 = (1,4,7,6,3)
class NucleotideNode{
public:
    NucleotideNode();
    ~NucleotideNode();
    void dfs();

public:
    int count;
    char base;
    NucleotideNode* children[8];
};

class NucleotideTree{
public:
    NucleotideTree(Options* opt);
    ~NucleotideTree();
    void addSeq(string seq);
    string getDominantPath(bool& reachedLeaf);

    static bool test();

private:
    Options* mOptions;
    NucleotideNode* mRoot;
};


#endif