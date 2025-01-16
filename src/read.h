#ifndef READ_H
#define READ_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include "sequence.h"
#include <vector>

using namespace std;

class Read{
public:
	Read(string* name, string* seq, string* strand, string* quality, bool phred64=false);
    Read(const char* name, const char* seq, const char* strand, const char* quality, bool phred64=false);
    ~Read();
	void print();
    void printFile(ofstream& file);
    Read* reverseComplement();
    string firstIndex();
    string lastIndex();
    // default is Q20
    int lowQualCount(int qual=20);
    int length();
    string toString();
    string toStringWithTag(string tag);
    void appendToString(string* target);
    void appendToStringWithTag(string* target, string tag);
    void resize(int len);
    void convertPhred64To33();
    void trimFront(int len);
    bool fixMGI();

public:
    static bool test();

private:


public:
	string* mName;
	string* mSeq;
	string* mStrand;
	string* mQuality;
};

class ReadPair{
public:
    ReadPair();
    ~ReadPair();
    void setPair(Read* left, Read* right);
    bool eof();

    // merge a pair, without consideration of seq error caused false INDEL
    Read* fastMerge();
public:
    Read* mLeft;
    Read* mRight;

public:
    static bool test();
};

struct ReadPack {
    Read** data;
    int count;
};

typedef struct ReadPack ReadPack;

#endif