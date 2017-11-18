#include "fastqreader.h"
#include "util.h"
#include <string.h>

FastqReader::FastqReader(string filename, bool hasQuality, bool phred64){
	mFilename = filename;
	mZipFile = NULL;
	mZipped = false;
	mPhred64 = phred64;
	mHasQuality = hasQuality;
	init();
}

FastqReader::~FastqReader(){
	close();
}

void FastqReader::init(){
	if (ends_with(mFilename, ".gz")){
		mZipFile = gzopen(mFilename.c_str(), "r");
		mZipped = true;
		Read* r = read();

		//test if it has quality line or not (fastq/fasta)
		if (r->mQuality[0] == '@')
			mHasQuality = false;
		delete r;
		gzrewind(mZipFile);
	}
	else if (isFastq(mFilename)){
		mFile.open(mFilename.c_str(), ifstream::in);
		mZipped = false;
	}
}

void FastqReader::getBytes(size_t& bytesRead, size_t& bytesTotal) {
	if(mZipped) {
		bytesRead = gzoffset(mZipFile);
	} else {
		bytesRead = mFile.tellg();
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
	const int maxLine = 1024;
	char line[maxLine];

	if(mZipped) {
		char* buf = NULL;
		memset(line, 0, maxLine);
		buf = gzgets(mZipFile, line, maxLine);

		// EOF or error, return an empty string
		if(!buf)
			return string();

		// optimize for short reads
		// reached end of line
		if(line[maxLine-2] == '\0' || line[maxLine-2] == '\n') {
			clearLineBreaks(line);
			return string(line);
		} else {
			string s(buf);
			while(buf) {
				memset(line, 0, maxLine);
				buf = gzgets(mZipFile, line, maxLine);
				//eof or error
				if(!buf)
					break;
				//reached end of line
				if(line[maxLine-2] == '\0' || line[maxLine-2] == '\n') {
					clearLineBreaks(buf);
					s.append(buf);
					break;
				} else {
					s.append(buf);
				}
			}
			return s;
		}
	}
	else {
		mFile.getline(line, maxLine);
		// optimize for short reads
		// reached end of line
		if(mFile.eof() || mFile.good()) {
			clearLineBreaks(line);
			return string(line);
		} else {
			string s(line);
			while(true) {
				if(mFile.eof())
					break;
				//clear fail bit
				mFile.clear();
				memset(line, 0, maxLine);
				mFile.getline(line, maxLine);
				if(mFile.eof() || mFile.good()) {
					clearLineBreaks(line);
					s.append(line);
					break;
				} else {
					// in case of some error happened, break it
					if(line[0] == '\0'){
						break;
					}
					s.append(line);
				}
			}
			return s;
		}
	}

	return string();
}

bool FastqReader::eof() {
	if (mZipped) {
		return gzeof(mZipFile);
	} else {
		return mFile.eof();
	}
}

Read* FastqReader::read(){
	if (mZipped){
		if (mZipFile == NULL)
			return NULL;
	}

	if(eof()) {
		return NULL;
	}

	string name = getLine();
	string sequence = getLine();
	string strand = getLine();

	if(name.empty() || sequence.empty() || strand.empty())
		return NULL;

	// WAR for FQ with no quality
	if (!mHasQuality){
		string quality = string(sequence.length(), 'K');
		return new Read(name, sequence, strand, quality, mPhred64);
	}
	else {
		string quality = getLine();
		if(quality.empty())
			return NULL;
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
		if (mFile.is_open()){
			mFile.close();
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

FastqReaderPair::FastqReaderPair(string leftName, string rightName, bool hasQuality, bool phred64){
	mLeft = new FastqReader(leftName, hasQuality, phred64);
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
	Read* r = mRight->read();
	if(!l || !r){
		return NULL;
	} else {
		return new ReadPair(l, r);
	}
}
