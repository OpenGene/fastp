
#include "fastareader.h"
#include "util.h"
#include <sstream>

FastaReader::FastaReader(string faFile, bool forceUpperCase)
{
    // Set locale and disable stdio synchronization to improve iostream performance
    // http://www.drdobbs.com/the-standard-librarian-iostreams-and-std/184401305
    // http://stackoverflow.com/questions/5166263/how-to-get-iostream-to-perform-better
    setlocale(LC_ALL,"C");
    ios_base::sync_with_stdio(false);

    mFastaFile = faFile;
    mForceUpperCase = forceUpperCase;
    if (is_directory(mFastaFile)) {
        string error_msg = "There is a problem with the provided fasta file: \'";
        error_msg.append(mFastaFile);
        error_msg.append("\' is a directory NOT a file...\n");
        throw invalid_argument(error_msg);
    }
    mFastaFileStream.open( mFastaFile.c_str(),ios::in);
    // verify that the file can be read
    if (!mFastaFileStream.is_open()) {
        string msg = "There is a problem with the provided fasta file: could NOT read ";
        msg.append(mFastaFile.c_str());
        msg.append("...\n");
        throw invalid_argument(msg);
    }

    char c;
    // seek to first contig
    while (mFastaFileStream.get(c) && c != '>') {
        if (mFastaFileStream.eof()) {
            break;
        }
    }
}

FastaReader::~FastaReader()
{
    if (mFastaFileStream.is_open()) {
        mFastaFileStream.close();
    }
}

void FastaReader::readNext()
{
    mCurrentID = "";
    mCurrentDescription = "";
    mCurrentSequence = "";
    bool foundHeader = false;
    
    char c;
    stringstream ssSeq;
    stringstream ssHeader;
    while(true){
        mFastaFileStream.get(c);
        if(c == '>' || mFastaFileStream.eof())
            break;
        else {
            if (foundHeader){
                if(mForceUpperCase && c>='a' && c<='z') {
                    c -= ('a' - 'A');
                }
                ssSeq << c;
            }
            else
                ssHeader << c;
        }

        string line = "";
        getline(mFastaFileStream,line,'\n');


        if(foundHeader == false) {
            ssHeader << line;
            foundHeader = true;
        }
        else {
            str_keep_valid_sequence(line, mForceUpperCase);
            ssSeq << line;
        }
    }
    mCurrentSequence = ssSeq.str();
    string header = ssHeader.str();

    mCurrentID = header;
}

bool FastaReader::hasNext() {
    return !mFastaFileStream.eof();
}

void FastaReader::readAll() {
    while(!mFastaFileStream.eof()){
        readNext();
        mAllContigs[mCurrentID] = mCurrentSequence;
    }
}

bool FastaReader::test(){
    FastaReader reader("testdata/tinyref.fa");
    reader.readAll();

    string contig1 = "GATCACAGGTCTATCACCCTATTAATTGGTATTTTCGTCTGGGGGGTGTGGAGCCGGAGCACCCTATGTCGCAGT";
    string contig2 = "GTCTGCACAGCCGCTTTCCACACAGAACCCCCCCCTCCCCCCGCTTCTGGCAAACCCCAAAAACAAAGAACCCTA";

    if(reader.mAllContigs.count("contig1") == 0 || reader.mAllContigs.count("contig2") == 0 )
        return false;

    if(reader.mAllContigs["contig1"] != contig1 || reader.mAllContigs["contig2"] != contig2 )
        return false;

    return true;

}



