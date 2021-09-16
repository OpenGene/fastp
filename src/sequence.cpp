#include "sequence.h"

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
    string str(origin->length(), 0);
    int len = origin->length();
    for(int c=0;c<origin->length();c++){
        char base = (*origin)[c];
        switch(base){
            case 'A':
            case 'a':
                str[len-c-1] = 'T';
                break;
            case 'T':
            case 't':
                str[len-c-1] = 'A';
                break;
            case 'C':
            case 'c':
                str[len-c-1] = 'G';
                break;
            case 'G':
            case 'g':
                str[len-c-1] = 'C';
                break;
            default:
                str[len-c-1] = 'N';
        }
    }
    return str;
}

Sequence Sequence::reverseComplement(){
    string*  str = new string(mStr->length(), 0);
    int len = mStr->length();
    for(int c=0;c<mStr->length();c++){
        char base = (*mStr)[c];
        switch(base){
            case 'A':
            case 'a':
                (*str)[len-c-1] = 'T';
                break;
            case 'T':
            case 't':
                (*str)[len-c-1] = 'A';
                break;
            case 'C':
            case 'c':
                (*str)[len-c-1] = 'G';
                break;
            case 'G':
            case 'g':
                (*str)[len-c-1] = 'C';
                break;
            default:
                (*str)[len-c-1] = 'N';
        }
    }
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