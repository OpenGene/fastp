#ifndef FASTA_READER_H
#define FASTA_READER_H

// includes
#include <cctype>
#include <clocale>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <map>

using namespace std;

class FastaReader
{
public:
    FastaReader(string fastaFile, bool forceUpperCase = true);
    ~FastaReader();
    bool hasNext();
    void readNext();
    void readAll();

    inline string currentID()
    {
        return mCurrentID;
    }

    inline string currentDescription()
    {
        return mCurrentDescription;
    }

    inline string currentSequence()
    {
        return mCurrentSequence;
    }

    inline map<string, string>& contigs() {
        return mAllContigs;
    }

    static bool test();


public:
    string mCurrentSequence;
    string mCurrentID ;
    string mCurrentDescription;
    map<string, string> mAllContigs;

private:
    bool readLine();
    bool endOfLine(char c);
    void setFastaSequenceIdDescription();

private:
    string mFastaFile;
    ifstream mFastaFileStream;
    bool mForceUpperCase;
};


#endif

