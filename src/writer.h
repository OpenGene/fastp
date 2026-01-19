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

#ifndef _WRITER_H
#define _WRITER_H

#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include <iostream>
#include <fstream>
#include "libdeflate.h"
#include "options.h"
#include <stdio.h>

using namespace std;

class Writer{
public:
	Writer(Options* opt, string filename, int compression, bool isSTDOUT = false);
	~Writer();
	bool isZipped();
	bool writeString(const string& str);
	bool writeString(string* str);
	bool write(const char* strdata, size_t size);
	void flush();
	string filename();

public:
	static bool test();

private:
	void init();
	void close();
	bool writeInternal(const char* strdata, size_t size);

private:
	string mFilename;
	libdeflate_compressor* mCompressor;
	//ofstream* mOutStream;
	FILE* mFP;
	bool mZipped;
	int mCompression;
	bool haveToClose;
	char* mBuffer;
	size_t mBufDataLen;
	size_t mBufSize;
	Options* mOptions;
	bool mSTDOUT;
};

#endif