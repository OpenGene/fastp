#ifndef FASTQ_READER_H
#define FASTQ_READER_H

#include <stdio.h>
#include <stdlib.h>
#include "read.h"
#ifdef DYNAMIC_ZLIB
  #include <zlib.h>
#else
  #include "zlib/zlib.h"
#endif
#include "common.h"
#include <iostream>
#include <fstream>

class FastqReader{
public:
	FastqReader(string filename, bool hasQuality = true, bool phred64=false);
	~FastqReader();
	bool isZipped();

	void getBytes(size_t& bytesRead, size_t& bytesTotal);

	//this function is not thread-safe
	//do not call read() of a same FastqReader object from different threads concurrently
	Read* read();
	bool eof();
	bool hasNoLineBreakAtEnd();

public:
	static bool isZipFastq(string filename);
	static bool isFastq(string filename);
	static bool test();

private:
	void init();
	void close();
	string getLine();
	void clearLineBreaks(char* line);
	void readToBuf();

private:
	string mFilename;
	gzFile mZipFile;
	FILE* mFile;
	bool mZipped;
	bool mHasQuality;
	bool mPhred64;
	char* mBuf;
	int mBufDataLen;
	int mBufUsedLen;
	bool mStdinMode;
	bool mHasNoLineBreakAtEnd;

};

class FastqReaderPair{
public:
	FastqReaderPair(FastqReader* left, FastqReader* right);
	FastqReaderPair(string leftName, string rightName, bool hasQuality = true, bool phred64 = false, bool interleaved = false);
	~FastqReaderPair();
	ReadPair* read();
public:
	FastqReader* mLeft;
	FastqReader* mRight;
	bool mInterleaved;
};

#endif
