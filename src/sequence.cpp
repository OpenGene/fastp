#include "sequence.h"
#include "simd.h"

Sequence::Sequence(){
}

Sequence::Sequence(string*  seq){
    mStr = seq;
}

Sequence::~Sequence(){
    if(mStr)
        delete mStr;
}

void Sequence::print(){
    std::cerr << *mStr;
}

int Sequence::length(){
    return mStr->length();
}

string Sequence::reverseComplement(string* origin) {
    int len = origin->length();
    string str(len, 0);
    fastp_simd::reverseComplement(origin->c_str(), &str[0], len);
    return str;
}

Sequence Sequence::reverseComplement(){
    int len = mStr->length();
    string* str = new string(len, 0);
    fastp_simd::reverseComplement(mStr->c_str(), &(*str)[0], len);
    return Sequence(str);
}

Sequence Sequence::operator~(){
    return reverseComplement();
}

bool Sequence::test(){
    Sequence s(new string("AAAATTTTCCCCGGGG"));
    Sequence rc = ~s;
    if (*(s.mStr) != "AAAATTTTCCCCGGGG" ){
        cerr << "Failed in reverseComplement() expect AAAATTTTCCCCGGGG, but get "<< *(s.mStr);
        return false;
    }
    if (*(rc.mStr) != "CCCCGGGGAAAATTTT" ){
        cerr << "Failed in reverseComplement() expect CCCCGGGGAAAATTTT, but get "<< *(rc.mStr);
        return false;
    }
    return true;
}