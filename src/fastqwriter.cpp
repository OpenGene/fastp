#include "fastqwriter.h"
#include "util.h"
#include <string.h>
#include "fastqreader.h"

FastqWriter::FastqWriter(string filename, int compression){
	mCompression = compression;
	mFilename = filename;
	mZipFile = NULL;
	mZipped = false;
	init();
}

FastqWriter::~FastqWriter(){
	close();
}

string FastqWriter::filename(){
	return mFilename;
}

void FastqWriter::init(){
	if (FastqReader::isZipFastq(mFilename)){
		mZipFile = gzopen(mFilename.c_str(), "w");
        gzsetparams(mZipFile, mCompression, Z_DEFAULT_STRATEGY);
		mZipped = true;
	}
	else {
		mFile.open(mFilename.c_str(), ifstream::out);
		mZipped = false;
	}
}

bool FastqWriter::writeLine(string linestr){
	const char* line = linestr.c_str();
	int size = linestr.length();
	int written;
	bool status;
	if(mZipped){
		written = gzwrite(mZipFile, line, size);
		gzputc(mZipFile, '\n');
		status = size == written;
	}
	else{
		mFile.write(line, size);
		mFile.put('\n');
		status = !mFile.fail();
	}

	return status;
}

bool FastqWriter::write(Read* r){
	bool success = true;
	success &= writeLine(r->mName);
	success &= writeLine(r->mSeq.mStr);
	success &= writeLine(r->mStrand);
	success &= writeLine(r->mQuality);
	return success;
}

void FastqWriter::close(){
	if (mZipped){
		if (mZipFile){
			gzflush(mZipFile, Z_FINISH);
			gzclose(mZipFile);
			mZipFile = NULL;
		}
	}
	else {
		if (mFile.is_open()){
			mFile.flush();
			mFile.close();
		}
	}
}

bool FastqWriter::isZipped(){
	return mZipped;
}

bool FastqWriter::test(){
	Read r("@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA",
		"CTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCATTCTATGATGTAGTTTTCCTTAGGAGGACATTTTTTACATGAAATTATTAACCTAAATAGAGTTGATC",
		"+",
		"AAAAA6EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEEEAEEEEEEEEEEEEEEEAEEE/EEEEEEEEEEAAEAEAAEEEAEEAA");
	FastqWriter writer("test.fq.gz");
	return writer.write(&r);
}

FastqWriterPair::FastqWriterPair(FastqWriter* left, FastqWriter* right){
	mLeft = left;
	mRight = right;
}

FastqWriterPair::FastqWriterPair(string leftName, string rightName){
	mLeft = new FastqWriter(leftName);
	mRight = new FastqWriter(rightName);
}

FastqWriterPair::~FastqWriterPair(){
	if(mLeft){
		delete mLeft;
		mLeft = NULL;
	}
	if(mRight){
		delete mRight;
		mRight = NULL;
	}
}

bool FastqWriterPair::write(ReadPair* pair){
	bool status = true;
	status &= mLeft->write(pair->mLeft);
	status &= mRight->write(pair->mRight);
	return status ;
}