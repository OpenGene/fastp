#include "read.h"
#include <sstream>
#include <cstring>
#include "util.h"

Read::Read(string* name, string* seq, string* strand, string* quality, bool phred64){
	mName = name;
	mSeq = seq;
	mStrand = strand;
	mQuality = quality;
	if(phred64)
		convertPhred64To33();
}

Read::Read(const char* name, const char* seq, const char* strand, const char* quality, bool phred64) {
	mName = new string(name);
	mSeq = new string(seq);
	mStrand = new string(strand);
	mQuality = new string(quality);
	if(phred64)
		convertPhred64To33();
}

Read::~Read(){
	if(mName)
		delete mName;
	if(mStrand)
		delete mStrand;
	if(mQuality)
		delete mQuality;
	if(mSeq)
		delete mSeq;
}

void Read::convertPhred64To33(){
	for(int i=0; i<mQuality->length(); i++) {
		(*mQuality)[i] = max(33, (*mQuality)[i] - (64-33));
	}
}

void Read::print(){
	std::cerr << *mName << endl;
	std::cerr << *(mSeq) << endl;
	std::cerr << *mStrand << endl;
	std::cerr << *mQuality << endl;
}

void Read::printFile(ofstream& file){
	file << *mName << endl;
	file << *mSeq << endl;
	file << *mStrand << endl;
	file << *mQuality << endl;
}

Read* Read::reverseComplement(){
	string seq = Sequence::reverseComplement(mSeq);
	string qual;
	qual.assign(mQuality->rbegin(), mQuality->rend());
	return new Read(mName->c_str(), seq.c_str(), "+", qual.c_str());
}

void Read::resize(int len) {
	if(len > length() || len<0)
		return ;
	mSeq->resize(len);
	mQuality->resize(len);
}
   
void Read::trimFront(int len){
	len = min(length()-1, len);
	mSeq->erase(0, len);
	mQuality->erase(0, len);
}

string Read::lastIndex(){
	int len = mName->length();
	if(len<5)
		return "";
	for(int i=len-3;i>=0;i--){
		if((*mName)[i]==':' || (*mName)[i]=='+'){
			return mName->substr(i+1, len-i);
		}
	}
	return "";
}

string Read::firstIndex(){
	int len = mName->length();
	int end = len;
	if(len<5)
		return "";
	for(int i=len-3;i>=0;i--){
		if((*mName)[i]=='+')
			end = i-1;
		if((*mName)[i]==':'){
			return mName->substr(i+1, end-i);
		}
	}
	return "";
}

int Read::lowQualCount(int qual){
	int count = 0;
	for(int q=0;q<mQuality->size();q++){
		if((*mQuality)[q] < qual + 33)
			count++;
	}
	return count;
}

int Read::length(){
	return mSeq->length();
}

string Read::toString() {
	return *mName + "\n" + *mSeq + "\n" + *mStrand + "\n" + *mQuality + "\n";
}

void Read::appendToString(string* target) {
	size_t size = mName->length() + mSeq->length() + mStrand->length() + mQuality->length() + 4;
	target->reserve(target->size() + size);
	target->append(*mName);
	target->push_back('\n');
	target->append(*mSeq);
	target->push_back('\n');
	target->append(*mStrand);
	target->push_back('\n');
	target->append(*mQuality);
	target->push_back('\n');
}

void Read::appendToStringWithTag(string* target, string tag) {
	size_t size = mName->length() + 1 + tag.length() + mSeq->length() + mStrand->length() + mQuality->length() + 4;
	target->reserve(target->size() + size);
	target->append(*mName);
	target->push_back(' ');
	target->append(tag);
	target->push_back('\n');
	target->append(*mSeq);
	target->push_back('\n');
	target->append(*mStrand);
	target->push_back('\n');
	target->append(*mQuality);
	target->push_back('\n');
}

string Read::toStringWithTag(string tag) {
	return *mName + " " + tag + "\n" + *mSeq + "\n" + *mStrand + "\n" + *mQuality + "\n";
}

bool Read::fixMGI() {
	int len = mName->length();
	if((*mName)[len-1]=='1' || (*mName)[len-1]=='2') {
		if((*mName)[len-2] == '/') {
			string* newName = new string(mName->substr(0, len-2) + " " + mName->substr(len-2, 2));
			delete mName;
			mName = newName;
			return true;
		}
	}
	return false;
}

bool Read::test(){
	Read r(new string("@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA"),
		new string("CTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCATTCTATGATGTAGTTTTCCTTAGGAGGACATTTTTTACATGAAATTATTAACCTAAATAGAGTTGATC"),
		new string("+"),
		new string("AAAAA6EEEEEEEEEEEEEEEEE#EEEEEEEEEEEEEEEEE/EEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEEEAEEEEEEEEEEEEEEEAEEE/EEEEEEEEEEAAEAEAAEEEAEEAA"));
	string idx = r.lastIndex();
	return idx == "GGTCCCGA";
}

ReadPair::ReadPair(){
	mLeft = NULL;
	mRight = NULL;
}

ReadPair::~ReadPair(){
	if(mLeft){
		delete mLeft;
		mLeft = NULL;
	}
	if(mRight){
		delete mRight;
		mRight = NULL;
	}
}

void ReadPair::setPair(Read* left, Read* right){
	mLeft = left;
	mRight = right;
}

bool ReadPair::eof(){
	return mLeft == NULL || mRight == NULL;
}

Read* ReadPair::fastMerge(){
	Read* rcRight = mRight->reverseComplement();
	int len1 = mLeft->length();
	int len2 = rcRight->length();
	// use the pointer directly for speed
	const char* str1 = mLeft->mSeq->c_str();
	const char* str2 = rcRight->mSeq->c_str();
	const char* qual1 = mLeft->mQuality->c_str();
	const char* qual2 = rcRight->mQuality->c_str();

	// we require at least 30 bp overlapping to merge a pair
	const int MIN_OVERLAP = 30;
	bool overlapped = false;
	int olen = MIN_OVERLAP;
	int diff = 0;
	// the diff count for 1 high qual + 1 low qual
	int lowQualDiff = 0;

	while(olen <= min(len1, len2)){
		diff = 0;
		lowQualDiff = 0;
		bool ok = true;
		int offset = len1 - olen;
		for(int i=0;i<olen;i++){
			if(str1[offset+i] != str2[i]){
				diff++;
				// one is >= Q30 and the other is <= Q15
				if((qual1[offset+i]>='?' && qual2[i]<='0') || (qual1[offset+i]<='0' && qual2[i]>='?')){
					lowQualDiff++;
				}
				// we disallow high quality diff, and only allow up to 3 low qual diff
				if(diff>lowQualDiff || lowQualDiff>=3){
					ok = false;
					break;
				}
			}
		}
		if(ok){
			overlapped = true;
			break;
		}
		olen++;
	}

	if(overlapped){
		int offset = len1 - olen;
		stringstream ss;
		ss << mLeft->mName << " merged offset:" << offset << " overlap:" << olen << " diff:" << diff;
		string mergedName = ss.str();
		string mergedSeq = mLeft->mSeq->substr(0, offset) + *(rcRight->mSeq);
		string mergedQual = mLeft->mQuality->substr(0, offset) + *(rcRight->mQuality);
		// quality adjuction and correction for low qual diff
		for(int i=0;i<olen;i++){
			if(str1[offset+i] != str2[i]){
				if(qual1[offset+i]>='?' && qual2[i]<='0'){
					mergedSeq[offset+i] = str1[offset+i];
					mergedQual[offset+i] = qual1[offset+i];
				} else {
					mergedSeq[offset+i] = str2[i];
					mergedQual[offset+i] = qual2[i];
				}
			} else {
				// add the quality of the pair to make a high qual
				mergedQual[offset+i] =  qual1[offset+i] + qual2[i] - 33;
			}
		}
		delete rcRight;
		return new Read(new string(mergedName), new string(mergedSeq), new string("+"), new string(mergedQual));
	}

	delete rcRight;
	return NULL;
}

bool ReadPair::test(){
	Read* left = new Read(new string("@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA"),
		new string("TTTTTTCTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCATTCTATGATGTAG"),
		new string("+"),
		new string("AAAAA6EEEEEEEEEEEEEEEEE#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"));
	Read* right = new Read(new string("@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA"),
		new string("AAAAAACTACACCATAGAATGACTATGAGTCTCATAAGAATGCACTCAACTAGTCATCACTCCTGTGTTTTCATAAGAAAAAACAGTGTTAGAGTCCAAGAG"),
		new string("+"),
		new string("AAAAA6EEEEE/EEEEEEEEEEE#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"));

	ReadPair pair;
	pair.setPair(left, right);
	Read* merged = pair.fastMerge();
	if(merged == NULL)
		return false;

	if(*(merged->mSeq) != "TTTTTTCTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCATTCTATGATGTAGTTTTTT")
		return false;

	return true;
}
