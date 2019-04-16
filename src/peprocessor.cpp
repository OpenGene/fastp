#include "peprocessor.h"
#include "fastqreader.h"
#include <iostream>
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"
#include "adaptertrimmer.h"
#include "basecorrector.h"
#include "jsonreporter.h"
#include "htmlreporter.h"
#include "polyx.h"

PairEndProcessor::PairEndProcessor(Options* opt){
    mOptions = opt;
    mProduceFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(opt);
    mOutStream1 = NULL;
    mZipFile1 = NULL;
    mOutStream2 = NULL;
    mZipFile2 = NULL;
    mUmiProcessor = new UmiProcessor(opt);

    int isizeBufLen = mOptions->insertSizeMax + 1;
    mInsertSizeHist = new long[isizeBufLen];
    memset(mInsertSizeHist, 0, sizeof(long)*isizeBufLen);
    mLeftWriter =  NULL;
    mRightWriter = NULL;
    mUnpairedLeftWriter =  NULL;
    mUnpairedRightWriter = NULL;
    mMergedWriter = NULL;
    mFailedWriter = NULL;

    mDuplicate = NULL;
    if(mOptions->duplicate.enabled) {
        mDuplicate = new Duplicate(mOptions);
    }
}

PairEndProcessor::~PairEndProcessor() {
    delete mInsertSizeHist;
    if(mDuplicate) {
        delete mDuplicate;
        mDuplicate = NULL;
    }
}

void PairEndProcessor::initOutput() {
    if(!mOptions->unpaired1.empty())
        mUnpairedLeftWriter = new WriterThread(mOptions, mOptions->unpaired1);

    if(!mOptions->unpaired2.empty() && mOptions->unpaired2 != mOptions->unpaired1)
        mUnpairedRightWriter = new WriterThread(mOptions, mOptions->unpaired2);

    if(mOptions->merge.enabled) {
        if(!mOptions->merge.out.empty())
            mMergedWriter = new WriterThread(mOptions, mOptions->merge.out);
    }

    if(!mOptions->failedOut.empty())
        mFailedWriter = new WriterThread(mOptions, mOptions->failedOut);

    if(mOptions->out1.empty())
        return;
    
    mLeftWriter = new WriterThread(mOptions, mOptions->out1);
    if(!mOptions->out2.empty())
        mRightWriter = new WriterThread(mOptions, mOptions->out2);
}

void PairEndProcessor::closeOutput() {
    if(mLeftWriter) {
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
    if(mRightWriter) {
        delete mRightWriter;
        mRightWriter = NULL;
    }
    if(mMergedWriter) {
        delete mMergedWriter;
        mMergedWriter = NULL;
    }
    if(mFailedWriter) {
        delete mFailedWriter;
        mFailedWriter = NULL;
    }
    if(mUnpairedLeftWriter) {
        delete mUnpairedLeftWriter;
        mLeftWriter = NULL;
    }
    if(mUnpairedRightWriter) {
        delete mUnpairedRightWriter;
        mRightWriter = NULL;
    }
}

void PairEndProcessor::initConfig(ThreadConfig* config) {
    if(mOptions->out1.empty())
        return;
    if(mOptions->split.enabled) {
        config->initWriterForSplit();
    }
}


bool PairEndProcessor::process(){
    if(!mOptions->split.enabled)
        initOutput();

    initPackRepository();
    std::thread producer(std::bind(&PairEndProcessor::producerTask, this));

    //TODO: get the correct cycles
    int cycle = 151;
    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        configs[t] = new ThreadConfig(mOptions, t, true);
        initConfig(configs[t]);
    }

    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&PairEndProcessor::consumerTask, this, configs[t]));
    }

    std::thread* leftWriterThread = NULL;
    std::thread* rightWriterThread = NULL;
    std::thread* unpairedLeftWriterThread = NULL;
    std::thread* unpairedRightWriterThread = NULL;
    std::thread* mergedWriterThread = NULL;
    std::thread* failedWriterThread = NULL;
    if(mLeftWriter)
        leftWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mLeftWriter));
    if(mRightWriter)
        rightWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mRightWriter));
    if(mUnpairedLeftWriter)
        unpairedLeftWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mUnpairedLeftWriter));
    if(mUnpairedRightWriter)
        unpairedRightWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mUnpairedRightWriter));
    if(mMergedWriter)
        mergedWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mMergedWriter));
    if(mFailedWriter)
        failedWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mFailedWriter));

    producer.join();
    for(int t=0; t<mOptions->thread; t++){
        threads[t]->join();
    }

    if(!mOptions->split.enabled) {
        if(leftWriterThread)
            leftWriterThread->join();
        if(rightWriterThread)
            rightWriterThread->join();
        if(unpairedLeftWriterThread)
            unpairedLeftWriterThread->join();
        if(unpairedRightWriterThread)
            unpairedRightWriterThread->join();
        if(mergedWriterThread)
            mergedWriterThread->join();
        if(failedWriterThread)
            failedWriterThread->join();
    }

    if(mOptions->verbose)
        loginfo("start to generate reports\n");

    // merge stats and filter results
    vector<Stats*> preStats1;
    vector<Stats*> postStats1;
    vector<Stats*> preStats2;
    vector<Stats*> postStats2;
    vector<FilterResult*> filterResults;
    for(int t=0; t<mOptions->thread; t++){
        preStats1.push_back(configs[t]->getPreStats1());
        postStats1.push_back(configs[t]->getPostStats1());
        preStats2.push_back(configs[t]->getPreStats2());
        postStats2.push_back(configs[t]->getPostStats2());
        filterResults.push_back(configs[t]->getFilterResult());
    }
    Stats* finalPreStats1 = Stats::merge(preStats1);
    Stats* finalPostStats1 = Stats::merge(postStats1);
    Stats* finalPreStats2 = Stats::merge(preStats2);
    Stats* finalPostStats2 = Stats::merge(postStats2);
    FilterResult* finalFilterResult = FilterResult::merge(filterResults);

    cerr << "Read1 before filtering:"<<endl;
    finalPreStats1->print();
    cerr << endl;
    cerr << "Read2 before filtering:"<<endl;
    finalPreStats2->print();
    cerr << endl;
    if(!mOptions->merge.enabled) {
        cerr << "Read1 after filtering:"<<endl;
        finalPostStats1->print();
        cerr << endl;
        cerr << "Read2 aftering filtering:"<<endl;
        finalPostStats2->print();
    } else {
        cerr << "Merged and filtered:"<<endl;
        finalPostStats1->print();
    }

    cerr << endl;
    cerr << "Filtering result:"<<endl;
    finalFilterResult->print();

    int* dupHist = NULL;
    double* dupMeanTlen = NULL;
    double* dupMeanGC = NULL;
    double dupRate = 0.0;
    if(mOptions->duplicate.enabled) {
        dupHist = new int[mOptions->duplicate.histSize];
        memset(dupHist, 0, sizeof(int) * mOptions->duplicate.histSize);
        dupMeanGC = new double[mOptions->duplicate.histSize];
        memset(dupMeanGC, 0, sizeof(double) * mOptions->duplicate.histSize);
        dupRate = mDuplicate->statAll(dupHist, dupMeanGC, mOptions->duplicate.histSize);
        cerr << endl;
        cerr << "Duplication rate: " << dupRate * 100.0 << "%" << endl;
    }

    // insert size distribution
    int peakInsertSize = getPeakInsertSize();
    cerr << endl;
    cerr << "Insert size peak (evaluated by paired-end reads): " << peakInsertSize << endl;

    if(mOptions->merge.enabled) {
        cerr << endl;
        cerr << "Read pairs merged: " << finalFilterResult->mMergedPairs << endl;
        if(finalPostStats1->getReads() > 0) {
            double postMergedPercent = 100.0 * finalFilterResult->mMergedPairs / finalPostStats1->getReads();
            double preMergedPercent = 100.0 * finalFilterResult->mMergedPairs / finalPreStats1->getReads();
            cerr << "% of original read pairs: " << preMergedPercent << "%" << endl;
            cerr << "% in reads after filtering: " << postMergedPercent << "%" << endl;
        }
        cerr << endl;
    }

    // make JSON report
    JsonReporter jr(mOptions);
    jr.setDupHist(dupHist, dupMeanGC, dupRate);
    jr.setInsertHist(mInsertSizeHist, peakInsertSize);
    jr.report(finalFilterResult, finalPreStats1, finalPostStats1, finalPreStats2, finalPostStats2);

    // make HTML report
    HtmlReporter hr(mOptions);
    hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.setInsertHist(mInsertSizeHist, peakInsertSize);
    hr.report(finalFilterResult, finalPreStats1, finalPostStats1, finalPreStats2, finalPostStats2);

    // clean up
    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }

    delete finalPreStats1;
    delete finalPostStats1;
    delete finalPreStats2;
    delete finalPostStats2;
    delete finalFilterResult;

    if(mOptions->duplicate.enabled) {
        delete[] dupHist;
        delete[] dupMeanGC;
    }

    delete[] threads;
    delete[] configs;

    if(leftWriterThread)
        delete leftWriterThread;
    if(rightWriterThread)
        delete rightWriterThread;
    if(unpairedLeftWriterThread)
        delete unpairedLeftWriterThread;
    if(unpairedRightWriterThread)
        delete unpairedRightWriterThread;
    if(mergedWriterThread)
        delete mergedWriterThread;
    if(failedWriterThread)
        delete failedWriterThread;

    if(!mOptions->split.enabled)
        closeOutput();

    return true;
}

int PairEndProcessor::getPeakInsertSize() {
    int peak = 0;
    long maxCount = -1;
    for(int i=0; i<mOptions->insertSizeMax; i++) {
        if(mInsertSizeHist[i] > maxCount) {
            peak = i;
            maxCount = mInsertSizeHist[i];
        }
    }
    return peak;
}

bool PairEndProcessor::processPairEnd(ReadPairPack* pack, ThreadConfig* config){
    string outstr1;
    string outstr2;
    string unpairedOut1;
    string unpairedOut2;
    string singleOutput;
    string mergedOutput;
    string failedOut;
    int readPassed = 0;
    int mergedCount = 0;
    for(int p=0;p<pack->count;p++){
        ReadPair* pair = pack->data[p];
        Read* or1 = pair->mLeft;
        Read* or2 = pair->mRight;

        int lowQualNum1 = 0;
        int nBaseNum1 = 0;
        int lowQualNum2 = 0;
        int nBaseNum2 = 0;

        // stats the original read before trimming
        config->getPreStats1()->statRead(or1);
        config->getPreStats2()->statRead(or2);

        // handling the duplication profiling
        if(mDuplicate)
            mDuplicate->statPair(or1, or2);

        // filter by index
        if(mOptions->indexFilter.enabled && mFilter->filterByIndex(or1, or2)) {
            delete pair;
            continue;
        }

        // umi processing
        if(mOptions->umi.enabled)
            mUmiProcessor->process(or1, or2);

        // trim in head and tail, and apply quality cut in sliding window
        int frontTrimmed1 = 0;
        int frontTrimmed2 = 0;
        Read* r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1, frontTrimmed1);
        Read* r2 = mFilter->trimAndCut(or2, mOptions->trim.front2, mOptions->trim.tail2, frontTrimmed2);

        if(r1 != NULL && r2!=NULL) {
            if(mOptions->polyGTrim.enabled)
                PolyX::trimPolyG(r1, r2, config->getFilterResult(), mOptions->polyGTrim.minLen);
        }
        bool isizeEvaluated = false;
        if(r1 != NULL && r2!=NULL && (mOptions->adapter.enabled || mOptions->correction.enabled)){
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire, mOptions->overlapDiffPercentLimit/100.0);
            // we only use thread 0 to evaluae ISIZE
            if(config->getThreadId() == 0) {
                statInsertSize(r1, r2, ov, frontTrimmed1, frontTrimmed2);
                isizeEvaluated = true;
            }
            if(mOptions->correction.enabled) {
                BaseCorrector::correctByOverlapAnalysis(r1, r2, config->getFilterResult(), ov);
            }
            if(mOptions->adapter.enabled) {
                bool trimmed = AdapterTrimmer::trimByOverlapAnalysis(r1, r2, config->getFilterResult(), ov, frontTrimmed1, frontTrimmed2);
                bool trimmed1 = trimmed;
                bool trimmed2 = trimmed;
                if(!trimmed){
                    if(mOptions->adapter.hasSeqR1)
                        trimmed1 = AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.sequence, false);
                    if(mOptions->adapter.hasSeqR2)
                        trimmed2 = AdapterTrimmer::trimBySequence(r2, config->getFilterResult(), mOptions->adapter.sequenceR2, true);
                }
                if(mOptions->adapter.hasFasta) {
                    AdapterTrimmer::trimByMultiSequences(r1, config->getFilterResult(), mOptions->adapter.seqsInFasta, false, !trimmed1);
                    AdapterTrimmer::trimByMultiSequences(r2, config->getFilterResult(), mOptions->adapter.seqsInFasta, true, !trimmed2);
                }
            }
        }

        if(config->getThreadId() == 0 && !isizeEvaluated && r1 != NULL && r2!=NULL) {
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire, mOptions->overlapDiffPercentLimit/100.0);
            statInsertSize(r1, r2, ov, frontTrimmed1, frontTrimmed2);
            isizeEvaluated = true;
        }

        if(r1 != NULL && r2!=NULL) {
            if(mOptions->polyXTrim.enabled)
                PolyX::trimPolyX(r1, r2, config->getFilterResult(), mOptions->polyXTrim.minLen);
        }

        if(r1 != NULL && r2!=NULL) {
            if( mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length())
                r1->resize(mOptions->trim.maxLen1);
            if( mOptions->trim.maxLen2 > 0 && mOptions->trim.maxLen2 < r2->length())
                r2->resize(mOptions->trim.maxLen2);
        }

        Read* merged = NULL;
        // merging mode
        bool mergeProcessed = false;
        if(mOptions->merge.enabled && r1 && r2) {
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire, mOptions->overlapDiffPercentLimit/100.0);
            if(ov.overlapped) {
                merged = OverlapAnalysis::merge(r1, r2, ov);
                int result = mFilter->passFilter(merged);
                config->addFilterResult(result, 2);
                if(result == PASS_FILTER) {
                    mergedOutput += merged->toString();
                    config->getPostStats1()->statRead(merged);
                    readPassed++;
                    mergedCount++;
                }
                delete merged;
                mergeProcessed = true;
            } else if(mOptions->merge.includeUnmerged){
                int result1 = mFilter->passFilter(r1);
                config->addFilterResult(result1, 1);
                if(result1 == PASS_FILTER) {
                    mergedOutput += r1->toString();
                    config->getPostStats1()->statRead(r1);
                }

                int result2 = mFilter->passFilter(r2);
                config->addFilterResult(result2, 1);
                if(result2 == PASS_FILTER) {
                    mergedOutput += r2->toString();
                    config->getPostStats1()->statRead(r2);
                }
                if(result1 == PASS_FILTER && result2 == PASS_FILTER )
                    readPassed++;
                mergeProcessed = true;
            }
        }

        if(!mergeProcessed) {

            int result1 = mFilter->passFilter(r1);
            int result2 = mFilter->passFilter(r2);

            config->addFilterResult(max(result1, result2), 2);

            if( r1 != NULL &&  result1 == PASS_FILTER && r2 != NULL && result2 == PASS_FILTER ) {
                
                if(mOptions->outputToSTDOUT && !mOptions->merge.enabled) {
                    singleOutput += r1->toString() + r2->toString();
                } else {
                    outstr1 += r1->toString();
                    outstr2 += r2->toString();
                }

                // stats the read after filtering
                if(!mOptions->merge.enabled) {
                    config->getPostStats1()->statRead(r1);
                    config->getPostStats2()->statRead(r2);
                }

                readPassed++;
            } else if( r1 != NULL &&  result1 == PASS_FILTER) {
                if(mUnpairedLeftWriter) {
                    unpairedOut1 += r1->toString();
                    if(mFailedWriter)
                        failedOut += or2->toStringWithTag(FAILED_TYPES[result2]);
                } else {
                    if(mFailedWriter) {
                        failedOut += or1->toStringWithTag("paired_read_is_failing");
                        failedOut += or2->toStringWithTag(FAILED_TYPES[result2]);
                    }
                }
            } else if( r2 != NULL && result2 == PASS_FILTER) {
                if(mUnpairedLeftWriter || mUnpairedRightWriter) {
                    unpairedOut2 += r2->toString();
                    if(mFailedWriter)
                        failedOut += or1->toStringWithTag(FAILED_TYPES[result1]);
                } else {
                    if(mFailedWriter) {
                        failedOut += or1->toStringWithTag(FAILED_TYPES[result1]);
                        failedOut += or2->toStringWithTag("paired_read_is_failing");
                    }
                }
            }
        }

        delete pair;
        // if no trimming applied, r1 should be identical to or1
        if(r1 != or1 && r1 != NULL)
            delete r1;
        // if no trimming applied, r1 should be identical to or1
        if(r2 != or2 && r2 != NULL)
            delete r2;
    }
    // if splitting output, then no lock is need since different threads write different files
    if(!mOptions->split.enabled) 
        mOutputMtx.lock();
    if(mOptions->outputToSTDOUT) {
        // STDOUT output
        // if it's merging mode, write the merged reads to STDOUT
        // otherwise write interleaved single output
        if(mOptions->merge.enabled)
            fwrite(mergedOutput.c_str(), 1, mergedOutput.length(), stdout);
        else
            fwrite(singleOutput.c_str(), 1, singleOutput.length(), stdout);
    } else if(mOptions->split.enabled) {
        // split output by each worker thread
        if(!mOptions->out1.empty())
            config->getWriter1()->writeString(outstr1);
        if(!mOptions->out2.empty())
            config->getWriter2()->writeString(outstr2);
    } 

    if(mMergedWriter && !mergedOutput.empty()) {
        // write merged data
        char* mdata = new char[mergedOutput.size()];
        memcpy(mdata, mergedOutput.c_str(), mergedOutput.size());
        mMergedWriter->input(mdata, mergedOutput.size());
    }

    if(mFailedWriter && !failedOut.empty()) {
        // write failed data
        char* fdata = new char[failedOut.size()];
        memcpy(fdata, failedOut.c_str(), failedOut.size());
        mFailedWriter->input(fdata, failedOut.size());
    }

    // normal output by left/right writer thread
    if(mRightWriter && mLeftWriter && (!outstr1.empty() || !outstr2.empty())) {
        // write PE
        char* ldata = new char[outstr1.size()];
        memcpy(ldata, outstr1.c_str(), outstr1.size());
        mLeftWriter->input(ldata, outstr1.size());

        char* rdata = new char[outstr2.size()];
        memcpy(rdata, outstr2.c_str(), outstr2.size());
        mRightWriter->input(rdata, outstr2.size());
    } else if(mLeftWriter && !singleOutput.empty()) {
        // write singleOutput
        char* ldata = new char[singleOutput.size()];
        memcpy(ldata, singleOutput.c_str(), singleOutput.size());
        mLeftWriter->input(ldata, singleOutput.size());
    }
    // output unpaired reads
    if (!unpairedOut1.empty() || !unpairedOut2.empty()) {
        if(mUnpairedLeftWriter && mUnpairedRightWriter) {
            // write PE
            char* unpairedData1 = new char[unpairedOut1.size()];
            memcpy(unpairedData1, unpairedOut1.c_str(), unpairedOut1.size());
            mUnpairedLeftWriter->input(unpairedData1, unpairedOut1.size());

            char* unpairedData2 = new char[unpairedOut2.size()];
            memcpy(unpairedData2, unpairedOut2.c_str(), unpairedOut2.size());
            mUnpairedRightWriter->input(unpairedData2, unpairedOut2.size());
        } else if(mUnpairedLeftWriter) {
            char* unpairedData = new char[unpairedOut1.size() + unpairedOut2.size() ];
            memcpy(unpairedData, unpairedOut1.c_str(), unpairedOut1.size());
            memcpy(unpairedData + unpairedOut1.size(), unpairedOut2.c_str(), unpairedOut2.size());
            mUnpairedLeftWriter->input(unpairedData, unpairedOut1.size() + unpairedOut2.size());
        }
    }

    if(!mOptions->split.enabled)
        mOutputMtx.unlock();

    if(mOptions->split.byFileLines)
        config->markProcessed(readPassed);
    else
        config->markProcessed(pack->count);

    if(mOptions->merge.enabled) {
        config->addMergedPairs(mergedCount);
    }

    delete pack->data;
    delete pack;

    return true;
}
    
void PairEndProcessor::statInsertSize(Read* r1, Read* r2, OverlapResult& ov, int frontTrimmed1, int frontTrimmed2) {
    int isize = mOptions->insertSizeMax;
    if(ov.overlapped) {
        if(ov.offset > 0)
            isize = r1->length() + r2->length() - ov.overlap_len + frontTrimmed1 + frontTrimmed2;
        else
            isize = ov.overlap_len + frontTrimmed1 + frontTrimmed2;
    }

    if(isize > mOptions->insertSizeMax)
        isize = mOptions->insertSizeMax;

    mInsertSizeHist[isize]++;
}

bool PairEndProcessor::processRead(Read* r, ReadPair* originalPair, bool reversed) {
    // do something here
    return true;
}

void PairEndProcessor::initPackRepository() {
    mRepo.packBuffer = new ReadPairPack*[PACK_NUM_LIMIT];
    memset(mRepo.packBuffer, 0, sizeof(ReadPairPack*)*PACK_NUM_LIMIT);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
    
}

void PairEndProcessor::destroyPackRepository() {
    delete mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void PairEndProcessor::producePack(ReadPairPack* pack){
    //std::unique_lock<std::mutex> lock(mRepo.mtx);
    /*while(((mRepo.writePos + 1) % PACK_NUM_LIMIT)
        == mRepo.readPos) {
        //mRepo.repoNotFull.wait(lock);
    }*/

    mRepo.packBuffer[mRepo.writePos] = pack;
    mRepo.writePos++;

    /*if (mRepo.writePos == PACK_NUM_LIMIT)
        mRepo.writePos = 0;*/

    //mRepo.repoNotEmpty.notify_all();
    //lock.unlock();
}

void PairEndProcessor::consumePack(ThreadConfig* config){
    ReadPairPack* data;
    //std::unique_lock<std::mutex> lock(mRepo.mtx);
    // buffer is empty, just wait here.
    /*while(mRepo.writePos % PACK_NUM_LIMIT == mRepo.readPos % PACK_NUM_LIMIT) {
        if(mProduceFinished){
            //lock.unlock();
            return;
        }
        //mRepo.repoNotEmpty.wait(lock);
    }*/

    mInputMtx.lock();
    while(mRepo.writePos <= mRepo.readPos) {
        usleep(1000);
        if(mProduceFinished) {
            mInputMtx.unlock();
            return;
        }
    }
    data = mRepo.packBuffer[mRepo.readPos];
    mRepo.readPos++;

    /*if (mRepo.readPos >= PACK_NUM_LIMIT)
        mRepo.readPos = 0;*/
    mInputMtx.unlock();
    //mRepo.readPos++;

    //lock.unlock();
    //mRepo.repoNotFull.notify_all();

    processPairEnd(data, config);

}

void PairEndProcessor::producerTask()
{
    if(mOptions->verbose)
        loginfo("start to load data");
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    ReadPair** data = new ReadPair*[PACK_SIZE];
    memset(data, 0, sizeof(ReadPair*)*PACK_SIZE);
    FastqReaderPair reader(mOptions->in1, mOptions->in2, true, mOptions->phred64, mOptions->interleavedInput);
    int count=0;
    bool needToBreak = false;
    while(true){
        ReadPair* read = reader.read();
        // TODO: put needToBreak here is just a WAR for resolve some unidentified dead lock issue 
        if(!read || needToBreak){
            // the last pack
            ReadPairPack* pack = new ReadPairPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            data = NULL;
            if(read) {
                delete read;
                read = NULL;
            }
            break;
        }
        data[count] = read;
        count++;
        // configured to process only first N reads
        if(mOptions->readsToProcess >0 && count + readNum >= mOptions->readsToProcess) {
            needToBreak = true;
        }
        if(mOptions->verbose && count + readNum >= lastReported + 1000000) {
            lastReported = count + readNum;
            string msg = "loaded " + to_string((lastReported/1000000)) + "M read pairs";
            loginfo(msg);
        }
        // a full pack
        if(count == PACK_SIZE || needToBreak){
            ReadPairPack* pack = new ReadPairPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            //re-initialize data for next pack
            data = new ReadPair*[PACK_SIZE];
            memset(data, 0, sizeof(ReadPair*)*PACK_SIZE);
            // if the consumer is far behind this producer, sleep and wait to limit memory usage
            while(mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT){
                slept++;
                usleep(1000);
            }
            readNum += count;
            // if the writer threads are far behind this producer, sleep and wait
            // check this only when necessary
            if(readNum % (PACK_SIZE * PACK_IN_MEM_LIMIT) == 0 && mLeftWriter) {
                while( (mLeftWriter && mLeftWriter->bufferLength() > PACK_IN_MEM_LIMIT) || (mRightWriter && mRightWriter->bufferLength() > PACK_IN_MEM_LIMIT) ){
                    slept++;
                    usleep(1000);
                }
            }
            // reset count to 0
            count = 0;
            // re-evaluate split size
            // TODO: following codes are commented since it may cause threading related conflicts in some systems
            /*if(mOptions->split.needEvaluation && !splitSizeReEvaluated && readNum >= mOptions->split.size) {
                splitSizeReEvaluated = true;
                // greater than the initial evaluation
                if(readNum >= 1024*1024) {
                    size_t bytesRead;
                    size_t bytesTotal;
                    reader.mLeft->getBytes(bytesRead, bytesTotal);
                    mOptions->split.size *=  (double)bytesTotal / ((double)bytesRead * (double) mOptions->split.number);
                    if(mOptions->split.size <= 0)
                        mOptions->split.size = 1;
                }
            }*/
        }
    }

    //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
    mProduceFinished = true;
    if(mOptions->verbose)
        loginfo("all reads loaded, start to monitor thread status");
    //lock.unlock();

    // if the last data initialized is not used, free it
    if(data != NULL)
        delete[] data;
}

void PairEndProcessor::consumerTask(ThreadConfig* config)
{
    while(true) {
        if(config->canBeStopped()){
            mFinishedThreads++;
            break;
        }
        while(mRepo.writePos <= mRepo.readPos) {
            if(mProduceFinished)
                break;
            usleep(1000);
        }
        //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
        if(mProduceFinished && mRepo.writePos == mRepo.readPos){
            mFinishedThreads++;
            if(mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
                loginfo(msg);
            }
            //lock.unlock();
            break;
        }
        if(mProduceFinished){
            if(mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " is processing the " + to_string(mRepo.readPos) + " / " + to_string(mRepo.writePos) + " pack";
                loginfo(msg);
            }
            consumePack(config);
            //lock.unlock();
        } else {
            //lock.unlock();
            consumePack(config);
        }
    }

    if(mFinishedThreads == mOptions->thread) {
        if(mLeftWriter)
            mLeftWriter->setInputCompleted();
        if(mRightWriter)
            mRightWriter->setInputCompleted();
        if(mUnpairedLeftWriter)
            mUnpairedLeftWriter->setInputCompleted();
        if(mUnpairedRightWriter)
            mUnpairedRightWriter->setInputCompleted();
        if(mMergedWriter)
            mMergedWriter->setInputCompleted();
        if(mFailedWriter)
            mFailedWriter->setInputCompleted();
    }
    
    if(mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
        loginfo(msg);
    }
}

void PairEndProcessor::writeTask(WriterThread* config)
{
    while(true) {
        if(config->isCompleted()){
            // last check for possible threading related issue
            config->output();
            break;
        }
        config->output();
    }

    if(mOptions->verbose) {
        string msg = config->getFilename() + " writer finished";
        loginfo(msg);
    }
}
