#include "threadconfig.h"
#include "util.h"

ThreadConfig::ThreadConfig(Options* opt, int threadId, bool paired){
    mOptions = opt;
    mThreadId = threadId;
    mWorkingSplit = threadId;
    mCurrentSplitReads = 0;
    mPreStats1 = new Stats(mOptions, false);
    mPostStats1 = new Stats(mOptions, false);
    if(paired){
        mPreStats2 = new Stats(mOptions, true);
        mPostStats2 = new Stats(mOptions, true);
    }
    else {
        mPreStats2 = NULL;
        mPostStats2 = NULL;
    }
    mWriter1 = NULL;
    mWriter2 = NULL;

    mFilterResult = new FilterResult(opt, paired);
    mCanBeStopped = false;
}

ThreadConfig::~ThreadConfig() {
    cleanup();
}

void ThreadConfig::cleanup() {
    if(mOptions->split.enabled && mOptions->split.byFileNumber)
        writeEmptyFilesForSplitting();
    deleteWriter();
}

void ThreadConfig::deleteWriter() {
    if(mWriter1 != NULL) {
        delete mWriter1;
        mWriter1 = NULL;
    }
    if(mWriter2 != NULL) {
        delete mWriter2;
        mWriter2 = NULL;
    }
}

void ThreadConfig::initWriter(string filename1) {
    deleteWriter();
    mWriter1 = new Writer(filename1, mOptions->compression);
}

void ThreadConfig::initWriter(string filename1, string filename2) {
    deleteWriter();
    mWriter1 = new Writer(filename1, mOptions->compression);
    mWriter2 = new Writer(filename2, mOptions->compression);
}

void ThreadConfig::initWriter(ofstream* stream) {
    deleteWriter();
    mWriter1 = new Writer(stream);
}

void ThreadConfig::initWriter(ofstream* stream1, ofstream* stream2) {
    deleteWriter();
    mWriter1 = new Writer(stream1);
    mWriter2 = new Writer(stream2);
}

void ThreadConfig::initWriter(gzFile gzfile) {
    deleteWriter();
    mWriter1 = new Writer(gzfile);
}

void ThreadConfig::initWriter(gzFile gzfile1, gzFile gzfile2) {
    deleteWriter();
    mWriter1 = new Writer(gzfile1);
    mWriter2 = new Writer(gzfile2);
}

void ThreadConfig::addFilterResult(int result, int readNum) {
    mFilterResult->addFilterResult(result, readNum);
}

void ThreadConfig::addMergedPairs(int pairs) {
    mFilterResult->addMergedPairs(pairs);
}

void ThreadConfig::initWriterForSplit() {
    if(mOptions->out1.empty())
        return ;

    // use 1-based naming
    string num = to_string(mWorkingSplit + 1);
    // padding for digits like 0001
    if(mOptions->split.digits > 0){
        while(num.size() < mOptions->split.digits)
            num = "0" + num;
    }

    string filename1 = joinpath(dirname(mOptions->out1), num + "." + basename(mOptions->out1));
    if(!mOptions->isPaired()) {
        initWriter(filename1);
    } else {
        string filename2 = joinpath(dirname(mOptions->out2), num + "." + basename(mOptions->out2));
        initWriter(filename1, filename2);
    }
}

void ThreadConfig::markProcessed(long readNum) {
    mCurrentSplitReads += readNum;
    if(!mOptions->split.enabled)
        return ;
    // if splitting is enabled, check whether current file is full
    if(mCurrentSplitReads >= mOptions->split.size) {
        // if it's splitting by file number, totally we cannot exceed split.number
        // if it's splitting by file lines, then we don't need to check
        if(mOptions->split.byFileLines || mWorkingSplit + mOptions->thread < mOptions->split.number ){
            mWorkingSplit += mOptions->thread;
            initWriterForSplit();
            mCurrentSplitReads = 0;
        } else {
            // this thread can be stoped now since all its tasks are done
            // only a part of threads have to deal with the remaining reads
            if(mOptions->split.number % mOptions->thread >0 
                && mThreadId >= mOptions->split.number % mOptions->thread)
                mCanBeStopped = true;
        }
    }
}

// if a task of writting N files is assigned to this thread, but the input file doesn't have so many reads to input
// write some empty files so it will not break following pipelines
void ThreadConfig::writeEmptyFilesForSplitting() {
    while(mWorkingSplit + mOptions->thread < mOptions->split.number) {
        mWorkingSplit += mOptions->thread;
            initWriterForSplit();
            mCurrentSplitReads = 0;
    }
}

bool ThreadConfig::canBeStopped() {
    return mCanBeStopped;
}