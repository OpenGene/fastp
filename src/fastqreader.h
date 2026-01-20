/*
MIT License

Copyright (c) 2021 Shifu Chen <chen@haplox.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef FASTQ_READER_H
#define FASTQ_READER_H

#include <stdio.h>
#include <stdlib.h>
#include "read.h"
#include "common.h"
#include <iostream>
#include <fstream>
#include "igzip_lib.h"
#include "readpool.h"

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
	void setReadPool(ReadPool* rp);

public:
	static bool isZipFastq(string filename);
	static bool isFastq(string filename);
	static bool test();

private:
	void init();
	void close();
	void getLine(string* line);
	void clearLineBreaks(char* line);
	void readToBuf();
	void readToBufIgzip();
	bool bufferFinished();

private:
	string mFilename;
	struct isal_gzip_header mGzipHeader;
	struct inflate_state mGzipState;
	unsigned char *mGzipInputBuffer;
	unsigned char *mGzipOutputBuffer;
	size_t mGzipInputBufferSize;
	size_t mGzipOutputBufferSize;
	size_t mGzipInputUsedBytes;
	FILE* mFile;
	bool mZipped;
	char* mFastqBuf;
	int mBufDataLen;
	int mBufUsedLen;
	bool mStdinMode;
	bool mHasNoLineBreakAtEnd;
	long mCounter;
	bool mHasQuality;
	bool mPhred64;
    ReadPool* mReadPool;

};

class FastqReaderPair{
public:
	FastqReaderPair(FastqReader* left, FastqReader* right);
	FastqReaderPair(string leftName, string rightName, bool hasQuality = true, bool phred64 = false, bool interleaved = false);
	~FastqReaderPair();
	void read(ReadPair* pair);
	bool eof();
public:
	FastqReader* mLeft;
	FastqReader* mRight;
	bool mInterleaved;
};

#endif
