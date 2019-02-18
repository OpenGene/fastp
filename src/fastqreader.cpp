#include "fastqreader.h"
#include "util.h"
#include <string.h>

#define FQ_BUF_SIZE (1<<20)

FastqReader::FastqReader(string filename, bool hasQuality, bool phred64){
	mFilename = filename;
	mZipFile = NULL;
	mZipped = false;
	mFile = NULL;
	mStdinMode = false;
	mPhred64 = phred64;
	mHasQuality = hasQuality;
	mBuf = new char[FQ_BUF_SIZE];
	mBufDataLen = 0;
	mBufUsedLen = 0;
	mHasNoLineBreakAtEnd = false;
	init();
}

FastqReader::~FastqReader(){
	close();
	delete mBuf;
}

bool FastqReader::hasNoLineBreakAtEnd() {
	return mHasNoLineBreakAtEnd;
}

void FastqReader::readToBuf() {
	if(mZipped) {
		mBufDataLen = gzread(mZipFile, mBuf, FQ_BUF_SIZE);
		if(mBufDataLen == -1) {
			cerr << "Error to read gzip file" << endl;
		}
	} else {
		mBufDataLen = fread(mBuf, 1, FQ_BUF_SIZE, mFile);
	}
	mBufUsedLen = 0;

	if(mBufDataLen < FQ_BUF_SIZE) {
		if(mBuf[mBufDataLen-1] != '\n')
			mHasNoLineBreakAtEnd = true;
	}
}

void FastqReader::init(){
	if (ends_with(mFilename, ".gz")){
		mZipFile = gzopen(mFilename.c_str(), "r");
		mZipped = true;
		gzrewind(mZipFile);
	}
	else {
		if(mFilename == "/dev/stdin") {
			mFile = stdin;
		}
		else
			mFile = fopen(mFilename.c_str(), "rb");
		if(mFile == NULL) {
			error_exit("Failed to open file: " + mFilename);
		}
		mZipped = false;
	}
	readToBuf();
}

void FastqReader::getBytes(size_t& bytesRead, size_t& bytesTotal) {
	if(mZipped) {
		bytesRead = gzoffset(mZipFile);
	} else {
		bytesRead = ftell(mFile);//mFile.tellg();
	}

	// use another ifstream to not affect current reader
	ifstream is(mFilename);
	is.seekg (0, is.end);
	bytesTotal = is.tellg();
}

void FastqReader::clearLineBreaks(char* line) {

	// trim \n, \r or \r\n in the tail
	int readed = strlen(line);
	if(readed >=2 ){
		if(line[readed-1] == '\n' || line[readed-1] == '\r'){
			line[readed-1] = '\0';
			if(line[readed-2] == '\r')
				line[readed-2] = '\0';
		}
	}
}

string FastqReader::getLine(){
	static int c=0;
	c++;
	int copied = 0;

	int start = mBufUsedLen;
	int end = start;

	while(end < mBufDataLen) {
		if(mBuf[end] != '\r' && mBuf[end] != '\n')
			end++;
		else
			break;
	}

	// this line well contained in this buf, or this is the last buf
	if(end < mBufDataLen || mBufDataLen < FQ_BUF_SIZE) {
		int len = end - start;
		string line(mBuf+start, len);

		// skip \n or \r
		end++;
		// handle \r\n
		if(end < mBufDataLen-1 && mBuf[end-1]=='\r' && mBuf[end] == '\n')
			end++;

		mBufUsedLen = end;

		return line;
	}

	// this line is not contained in this buf, we need to read new buf
	string str(mBuf+start, mBufDataLen - start);

	while(true) {
		readToBuf();
		start = 0;
		end = 0;
		while(end < mBufDataLen) {
			if(mBuf[end] != '\r' && mBuf[end] != '\n')
				end++;
			else
				break;
		}
		// this line well contained in this buf, we need to read new buf
		if(end < mBufDataLen || mBufDataLen < FQ_BUF_SIZE) {
			int len = end - start;
			str.append(mBuf+start, len);

			// skip \n or \r
			end++;
			// handle \r\n
			if(end < mBufDataLen-1 && mBuf[end] == '\n')
				end++;

			mBufUsedLen = end;
			return str;
		}
		// even this new buf is not enough, although impossible
		str.append(mBuf+start, mBufDataLen);
	}

	return string();
}

bool FastqReader::eof() {
	if (mZipped) {
		return gzeof(mZipFile);
	} else {
		return feof(mFile);//mFile.eof();
	}
}

Read* FastqReader::read(){
	if (mZipped){
		if (mZipFile == NULL)
			return NULL;
	}

	if(mBufUsedLen >= mBufDataLen && eof()) {
		return NULL;
	}

	string name = getLine();
	// name should start with @
	while((name.empty() && !(mBufUsedLen >= mBufDataLen && eof())) || (!name.empty() && name[0]!='@')){
		name = getLine();
	}

	if(name.empty())
		return NULL;

	string sequence = getLine();
	string strand = getLine();

	// WAR for FQ with no quality
	if (!mHasQuality){
		string quality = string(sequence.length(), 'K');
		return new Read(name, sequence, strand, quality, mPhred64);
	}
	else {
		string quality = getLine();
		if(quality.length() != sequence.length()) {
			cerr << "ERROR: sequence and quality have different length:" << endl;
			cerr << name << endl;
			cerr << sequence << endl;
			cerr << strand << endl;
			cerr << quality << endl;
			return NULL;
		}
		return new Read(name, sequence, strand, quality, mPhred64);
	}

	return NULL;
}

void FastqReader::close(){
	if (mZipped){
		if (mZipFile){
			gzclose(mZipFile);
			mZipFile = NULL;
		}
	}
	else {
		if (mFile){
			fclose(mFile);//mFile.close();
			mFile = NULL;
		}
	}
}

bool FastqReader::isZipFastq(string filename) {
	if (ends_with(filename, ".fastq.gz"))
		return true;
	else if (ends_with(filename, ".fq.gz"))
		return true;
	else if (ends_with(filename, ".fasta.gz"))
		return true;
	else if (ends_with(filename, ".fa.gz"))
		return true;
	else
		return false;
}

bool FastqReader::isFastq(string filename) {
	if (ends_with(filename, ".fastq"))
		return true;
	else if (ends_with(filename, ".fq"))
		return true;
	else if (ends_with(filename, ".fasta"))
		return true;
	else if (ends_with(filename, ".fa"))
		return true;
	else
		return false;
}

bool FastqReader::isZipped(){
	return mZipped;
}

bool FastqReader::test(){
	FastqReader reader1("testdata/R1.fq");
	FastqReader reader2("testdata/R1.fq.gz");
	Read* r1 = NULL;
	Read* r2 = NULL;
	while(true){
		r1=reader1.read();
		r2=reader2.read();
		if(r1 == NULL || r2 == NULL)
			break;
		if(r1->mSeq.mStr != r2->mSeq.mStr){
			return false;
		}
		delete r1;
		delete r2;
	}
	return true;
}

FastqReaderPair::FastqReaderPair(FastqReader* left, FastqReader* right){
	mLeft = left;
	mRight = right;
}

FastqReaderPair::FastqReaderPair(string leftName, string rightName, bool hasQuality, bool phred64, bool interleaved){
	mInterleaved = interleaved;
	mLeft = new FastqReader(leftName, hasQuality, phred64);
	if(mInterleaved)
		mRight = NULL;
	else
		mRight = new FastqReader(rightName, hasQuality, phred64);
}

FastqReaderPair::~FastqReaderPair(){
	if(mLeft){
		delete mLeft;
		mLeft = NULL;
	}
	if(mRight){
		delete mRight;
		mRight = NULL;
	}
}

ReadPair* FastqReaderPair::read(){
	Read* l = mLeft->read();
	Read* r = NULL;
	if(mInterleaved)
		r = mLeft->read();
	else
		r = mRight->read();
	if(!l || !r){
		return NULL;
	} else {
		return new ReadPair(l, r);
	}
}
