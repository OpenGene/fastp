#ifndef _WRITER_H
#define _WRITER_H

#include <stdio.h>
#include <stdlib.h>
#ifdef DYNAMIC_ZLIB
  #include <zlib.h>
#else
  #include "zlib/zlib.h"
#endif
#include "common.h"
#include <iostream>
#include <fstream>

using namespace std;

class Writer{
public:
	Writer(string filename, int compression = 3);
	Writer(ofstream* stream);
	Writer(gzFile gzfile);
	~Writer();
	bool isZipped();
	bool writeString(string& s);
	bool writeLine(string& linestr);
	bool write(char* strdata, size_t size);
	string filename();

public:
	static bool test();

private:
	void init();
	void close();

private:
	string mFilename;
	gzFile mZipFile;
	ofstream* mOutStream;
	bool mZipped;
	int mCompression;
	bool haveToClose;
};

#endif
