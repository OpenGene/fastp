#include <stdlib.h>
#include "filterresult.h"

FilterResult::FilterResult(Options* opt, bool paired){
    mOptions = opt;
    mPaired = paired;
    mTrimmedAdapterRead = 0;
    mTrimmedAdapterBases = 0;
    for(int i=0; i<FILTER_RESULT_TYPES; i++) {
        mFilterReadStats[i] = 0;
    }
}

FilterResult::~FilterResult() {
}

void FilterResult::addFilterResult(int result) {
    if(result < PASS_FILTER || result >= FILTER_RESULT_TYPES)
        return ;
    mFilterReadStats[result]++;
}

FilterResult* FilterResult::merge(vector<FilterResult*>& list) {
    if(list.size() == 0)
        return NULL;
    FilterResult* result = new FilterResult(list[0]->mOptions, list[0]->mPaired);

    long* target = result->getFilterReadStats();
    
    for(int i=0; i<list.size(); i++) {
        long* current = list[i]->getFilterReadStats();
        for(int j=0; j<FILTER_RESULT_TYPES; j++) {
            target[j] += current[j];
        }
        result->mTrimmedAdapterRead += list[i]->mTrimmedAdapterRead;
        result->mTrimmedAdapterBases += list[i]->mTrimmedAdapterBases;
    }
    return result;
}

void FilterResult::addAdapterTrimmed(string adapter) {
    mTrimmedAdapterRead++;
    mTrimmedAdapterBases += adapter.length();
}

void FilterResult::addAdapterTrimmed(string adapter1, string adapter2) {
    mTrimmedAdapterRead++;
    mTrimmedAdapterBases += adapter1.length() + adapter2.length();
}

void FilterResult::print() {
    cout << (mPaired?"Read pairs":"Reads") << " passed filter: " << mFilterReadStats[PASS_FILTER] << endl;
    cout << (mPaired?"Read pairs":"Reads") << " failed due to low quality: " << mFilterReadStats[FAIL_QUALITY] << endl;
    cout << (mPaired?"Read pairs":"Reads") << " failed due to too many N: " << mFilterReadStats[FAIL_N_BASE] << endl;
    cout << (mPaired?"Read pairs":"Reads") << " failed due to too short: " << mFilterReadStats[FAIL_LENGTH] << endl;
    cout << (mPaired?"Read pairs":"Reads") << " with adapter trimmed: " << mTrimmedAdapterRead << endl;
    cout << "Bases" << " trimmed for adapters: " << mTrimmedAdapterBases << endl;
}