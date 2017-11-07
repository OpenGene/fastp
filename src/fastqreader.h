#ifndef FASTQ_READER_H
#define FASTQ_READER_H

#include <stdio.h>
#include <stdlib.h>
#include "read.h"
#include "zlib/zlib.h"
#include "common.h"
#include <iostream>
#include <fstream>

class FastqReader{
public:
	FastqReader(string filename, bool hasQuality = true);
	~FastqReader();
	bool isZipped();

	void getBytes(size_t& bytesRead, size_t& bytesTotal);

	//this function is not thread-safe
	//do not call read() of a same FastqReader object from different threads concurrently
	Read* read();

public:
	static bool isZipFastq(string filename);
	static bool isFastq(string filename);
	static bool test();

private:
	void init();
	void close();
	bool getLine(char* line, int maxLine);

private:
	string mFilename;
	gzFile mZipFile;
	ifstream mFile;
	bool mZipped;
	bool mHasQuality;

};

class FastqReaderPair{
public:
	FastqReaderPair(FastqReader* left, FastqReader* right);
	FastqReaderPair(string leftName, string rightName);
	~FastqReaderPair();
	ReadPair* read();
public:
	FastqReader* mLeft;
	FastqReader* mRight;
};

#endif