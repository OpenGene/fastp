#ifndef FASTQ_WRITER_H
#define FASTQ_WRITER_H

#include <stdio.h>
#include <stdlib.h>
#include "read.h"
#include "zlib/zlib.h"
#include "common.h"
#include <iostream>
#include <fstream>

class FastqWriter{
public:
	FastqWriter(string filename, int compression = 3);
	~FastqWriter();
	bool isZipped();

	//this function is not thread-safe
	bool write(Read* r);
	void close();

	string filename();

public:
	static bool test();

private:
	void init();
	bool writeLine(string linestr);

private:
	string mFilename;
	gzFile mZipFile;
	ofstream mFile;
	bool mZipped;
	int mCompression;

};

class FastqWriterPair{
public:
	FastqWriterPair(FastqWriter* left, FastqWriter* right);
	FastqWriterPair(string leftName, string rightName);
	~FastqWriterPair();
	bool write(ReadPair* pair);
public:
	FastqWriter* mLeft;
	FastqWriter* mRight;
};

#endif