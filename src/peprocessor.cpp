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
    mLeftReaderFinished = false;
    mRightReaderFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(opt);
    mUmiProcessor = new UmiProcessor(opt);

    int isizeBufLen = mOptions->insertSizeMax + 1;
    mInsertSizeHist = new atomic_long[isizeBufLen];
    memset(mInsertSizeHist, 0, sizeof(atomic_long)*isizeBufLen);
    mLeftWriter =  NULL;
    mRightWriter = NULL;
    mUnpairedLeftWriter =  NULL;
    mUnpairedRightWriter = NULL;
    mMergedWriter = NULL;
    mFailedWriter = NULL;
    mOverlappedWriter = NULL;
    shouldStopReading = false;

    mDuplicate = NULL;
    if(mOptions->duplicate.enabled) {
        mDuplicate = new Duplicate(mOptions);
    }

    mLeftPackReadCounter = 0;
    mRightPackReadCounter = 0;
    mPackProcessedCounter = 0;

    mLeftReadPool = new ReadPool(mOptions);
    mRightReadPool = new ReadPool(mOptions);
}

PairEndProcessor::~PairEndProcessor() {
    delete mInsertSizeHist;
    if(mDuplicate) {
        delete mDuplicate;
        mDuplicate = NULL;
    }
    if(mLeftReadPool) {
        delete mLeftReadPool;
        mLeftReadPool = NULL;
    }
    if(mRightReadPool) {
        delete mRightReadPool;
        mRightReadPool = NULL;
    }
    delete[] mLeftInputLists;
    delete[] mRightInputLists;
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

    if(!mOptions->overlappedOut.empty())
        mOverlappedWriter = new WriterThread(mOptions, mOptions->overlappedOut);

    if(mOptions->out1.empty() && !mOptions->outputToSTDOUT)
        return;
    
    mLeftWriter = new WriterThread(mOptions, mOptions->out1, mOptions->outputToSTDOUT);
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
    if(mOverlappedWriter) {
        delete mOverlappedWriter;
        mOverlappedWriter = NULL;
    }
    if(mUnpairedLeftWriter) {
        delete mUnpairedLeftWriter;
        mUnpairedLeftWriter = NULL;
    }
    if(mUnpairedRightWriter) {
        delete mUnpairedRightWriter;
        mUnpairedRightWriter = NULL;
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

    std::thread* readerLeft = NULL;
    std::thread* readerRight = NULL;
    std::thread* readerInterveleaved = NULL;

    mLeftInputLists = new SingleProducerSingleConsumerList<ReadPack*>*[mOptions->thread];
    mRightInputLists = new SingleProducerSingleConsumerList<ReadPack*>*[mOptions->thread];

    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        mLeftInputLists[t] = new SingleProducerSingleConsumerList<ReadPack*>();
        mRightInputLists[t] = new SingleProducerSingleConsumerList<ReadPack*>();
        configs[t] = new ThreadConfig(mOptions, t, true);
        configs[t]->setInputListPair(mLeftInputLists[t], mRightInputLists[t]);
        initConfig(configs[t]);
    }

    if(mOptions->interleavedInput)
        readerInterveleaved= new std::thread(std::bind(&PairEndProcessor::interleavedReaderTask, this));
    else {
        readerLeft = new std::thread(std::bind(&PairEndProcessor::readerTask, this, true));
        readerRight = new std::thread(std::bind(&PairEndProcessor::readerTask, this, false));
    }

    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&PairEndProcessor::processorTask, this, configs[t]));
    }

    std::thread* leftWriterThread = NULL;
    std::thread* rightWriterThread = NULL;
    std::thread* unpairedLeftWriterThread = NULL;
    std::thread* unpairedRightWriterThread = NULL;
    std::thread* mergedWriterThread = NULL;
    std::thread* failedWriterThread = NULL;
    std::thread* overlappedWriterThread = NULL;
    if(mLeftWriter)
        leftWriterThread = new std::thread(std::bind(&PairEndProcessor::writerTask, this, mLeftWriter));
    if(mRightWriter)
        rightWriterThread = new std::thread(std::bind(&PairEndProcessor::writerTask, this, mRightWriter));
    if(mUnpairedLeftWriter)
        unpairedLeftWriterThread = new std::thread(std::bind(&PairEndProcessor::writerTask, this, mUnpairedLeftWriter));
    if(mUnpairedRightWriter)
        unpairedRightWriterThread = new std::thread(std::bind(&PairEndProcessor::writerTask, this, mUnpairedRightWriter));
    if(mMergedWriter)
        mergedWriterThread = new std::thread(std::bind(&PairEndProcessor::writerTask, this, mMergedWriter));
    if(mFailedWriter)
        failedWriterThread = new std::thread(std::bind(&PairEndProcessor::writerTask, this, mFailedWriter));
    if(mOverlappedWriter)
        overlappedWriterThread = new std::thread(std::bind(&PairEndProcessor::writerTask, this, mOverlappedWriter));

    if(readerInterveleaved) {
        readerInterveleaved->join();
    } else {
        readerLeft->join();
        readerRight->join();
    }
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
        if(overlappedWriterThread)
            overlappedWriterThread->join();
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
        cerr << "Read2 after filtering:"<<endl;
        finalPostStats2->print();
    } else {
        cerr << "Merged and filtered:"<<endl;
        finalPostStats1->print();
    }

    cerr << endl;
    cerr << "Filtering result:"<<endl;
    finalFilterResult->print();

    double dupRate = 0.0;
    if(mOptions->duplicate.enabled) {
        dupRate = mDuplicate->getDupRate();
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
    jr.setDup(dupRate);
    jr.setInsertHist(mInsertSizeHist, peakInsertSize);
    jr.report(finalFilterResult, finalPreStats1, finalPostStats1, finalPreStats2, finalPostStats2);

    // make HTML report
    HtmlReporter hr(mOptions);
    hr.setDup(dupRate);
    hr.setInsertHist(mInsertSizeHist, peakInsertSize);
    hr.report(finalFilterResult, finalPreStats1, finalPostStats1, finalPreStats2, finalPostStats2);

    // clean up
    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }

    if(readerInterveleaved) {
        delete readerInterveleaved;
    } else {
        delete readerLeft;
        delete readerRight;
    }

    delete finalPreStats1;
    delete finalPostStats1;
    delete finalPreStats2;
    delete finalPostStats2;
    delete finalFilterResult;

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
    if(overlappedWriterThread)
        delete overlappedWriterThread;

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

void PairEndProcessor::recycleToPool1(int tid, Read* r) {
    // failed to recycle, then delete it
    if(!mLeftReadPool->input(tid, r))
        delete r;
}

void PairEndProcessor::recycleToPool2(int tid, Read* r) {
    // failed to recycle, then delete it
    if(!mRightReadPool->input(tid, r))
        delete r;
}

bool PairEndProcessor::processPairEnd(ReadPack* leftPack, ReadPack* rightPack, ThreadConfig* config){
    if(leftPack->count != rightPack->count) {
        cerr << endl;
        cerr << "WARNING: different read numbers of the " << mPackProcessedCounter << " pack" << endl;
        cerr << "Read1 pack size: " << leftPack->count << endl;
        cerr << "Read2 pack size: " << rightPack->count << endl;
        cerr << "Ignore the unmatched reads" << endl << endl;
        shouldStopReading = true;
    }
    int tid = config->getThreadId();

    // build output on stack strings, move to heap only when handing off to writers
    string outstr1, outstr2, unpairedOut1, unpairedOut2;
    string singleOutput, mergedOutput, failedOut, overlappedOut;
    // reserve capacity for main outputs to avoid repeated reallocation
    const size_t estimatedCapacity = leftPack->count * 320;
    outstr1.reserve(estimatedCapacity);
    outstr2.reserve(estimatedCapacity);

    int readPassed = 0;
    int mergedCount = 0;
    for(int p=0;p<leftPack->count && p<rightPack->count;p++){
        Read* or1 = leftPack->data[p];
        Read* or2 = rightPack->data[p];

        int lowQualNum1 = 0;
        int nBaseNum1 = 0;
        int lowQualNum2 = 0;
        int nBaseNum2 = 0;

        // stats the original read before trimming
        config->getPreStats1()->statRead(or1);
        config->getPreStats2()->statRead(or2);

        // handling the duplication profiling
        bool dedupOut = false;
        if(mDuplicate) {
            bool isDup = mDuplicate->checkPair(or1, or2);
            if(mOptions->duplicate.dedup && isDup)
                dedupOut = true;
        }

        // filter by index
        if(mOptions->indexFilter.enabled && mFilter->filterByIndex(or1, or2)) {
            recycleToPool1(tid, or1);
            or1 = NULL;
            recycleToPool2(tid, or2);
            or2 = NULL;
            continue;
        }

        // fix MGI
        if(mOptions->fixMGI) {
            or1->fixMGI();
            or2->fixMGI();
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
        bool isAdapterDimer = false;
        if(r1 != NULL && r2!=NULL && (mOptions->adapter.enabled || mOptions->correction.enabled)){
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire, mOptions->overlapDiffPercentLimit/100.0, mOptions->adapter.allowGapOverlapTrimming);
            // we only use thread 0 to evaluae ISIZE
            if(config->getThreadId() == 0) {
                statInsertSize(r1, r2, ov, frontTrimmed1, frontTrimmed2);
                isizeEvaluated = true;
            }
            if(mOptions->correction.enabled && !ov.hasGap) {
                // no gap allowed for overlap correction
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
                    trimmed1 |= AdapterTrimmer::trimByMultiSequences(r1, config->getFilterResult(), mOptions->adapter.seqsInFasta, false, !trimmed1);
                    trimmed2 |= AdapterTrimmer::trimByMultiSequences(r2, config->getFilterResult(), mOptions->adapter.seqsInFasta, true, !trimmed2);
                }

                // Check for adapter dimer: both reads shorter than threshold after adapter trimming
                // AND adapters were detected in at least one of the reads (requires evidence)
                if(r1 != NULL && r2 != NULL && (trimmed1 || trimmed2) &&
                   r1->length() <= mOptions->adapter.dimerMaxLen &&
                   r2->length() <= mOptions->adapter.dimerMaxLen) {
                    isAdapterDimer = true;
                }
            }
        }

        if(r1 != NULL && r2!=NULL && mOverlappedWriter) {
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire, 0);
            if(ov.overlapped) {
                Read* overlappedRead = new Read(new string(*r1->mName), new string(r1->mSeq->substr(max(0,ov.offset)), ov.overlap_len), new string(*r1->mStrand), new string(r1->mQuality->substr(max(0,ov.offset)), ov.overlap_len));
                overlappedRead->appendToString(&overlappedOut);
                recycleToPool1(tid, overlappedRead);
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
                    merged->appendToString(&mergedOutput);
                    config->getPostStats1()->statRead(merged);
                    readPassed++;
                    mergedCount++;
                }
                recycleToPool1(tid, merged);
                mergeProcessed = true;
            } else if(mOptions->merge.includeUnmerged){
                int result1 = mFilter->passFilter(r1);
                int result2 = mFilter->passFilter(r2);

                if(isAdapterDimer) {
                    result1 = FAIL_ADAPTER_DIMER;
                    result2 = FAIL_ADAPTER_DIMER;
                }

                config->addFilterResult(result1, 1);
                if(result1 == PASS_FILTER && !dedupOut) {
                    r1->appendToString(&mergedOutput);
                    config->getPostStats1()->statRead(r1);
                }

                config->addFilterResult(result2, 1);
                if(result2 == PASS_FILTER && !dedupOut) {
                    r2->appendToString(&mergedOutput);
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

            if(isAdapterDimer) {
                result1 = FAIL_ADAPTER_DIMER;
                result2 = FAIL_ADAPTER_DIMER;
            }

            config->addFilterResult(max(result1, result2), 2);

            if(!dedupOut) {

                if( r1 != NULL &&  result1 == PASS_FILTER && r2 != NULL && result2 == PASS_FILTER ) {

                    if(mOptions->outputToSTDOUT && !mOptions->merge.enabled) {
                        r1->appendToString(&singleOutput);
                        r2->appendToString(&singleOutput);
                    } else {
                        r1->appendToString(&outstr1);
                        r2->appendToString(&outstr2);
                    }

                    // stats the read after filtering
                    if(!mOptions->merge.enabled) {
                        config->getPostStats1()->statRead(r1);
                        config->getPostStats2()->statRead(r2);
                    }

                    readPassed++;
                } else if( r1 != NULL &&  result1 == PASS_FILTER) {
                    if(mUnpairedLeftWriter) {
                        r1->appendToString(&unpairedOut1);
                        if(mFailedWriter)
                            or2->appendToStringWithTag(&failedOut, FAILED_TYPES[result2]);
                    } else {
                        if(mFailedWriter) {
                            or1->appendToStringWithTag(&failedOut, "paired_read_is_failing");
                            or2->appendToStringWithTag(&failedOut, FAILED_TYPES[result2]);
                        }
                    }
                } else if( r2 != NULL && result2 == PASS_FILTER) {
                    if(mUnpairedRightWriter) {
                        r2->appendToString(&unpairedOut2);
                        if(mFailedWriter)
                            or1->appendToStringWithTag(&failedOut,FAILED_TYPES[result1]);
                    } else if(mUnpairedLeftWriter) {
                        r2->appendToString(&unpairedOut1);
                        if(mFailedWriter)
                            or1->appendToStringWithTag(&failedOut,FAILED_TYPES[result1]);
                    }  else {
                        if(mFailedWriter) {
                            or1->appendToStringWithTag(&failedOut, FAILED_TYPES[result1]);
                            or2->appendToStringWithTag(&failedOut, "paired_read_is_failing");
                        }
                    }
                }
            }
        }

        // if no trimming applied, r1 should be identical to or1
        if(r1 != or1 && r1 != NULL) {
            recycleToPool1(tid, r1);
            r1 = NULL;
        }
        // if no trimming applied, r1 should be identical to or1
        if(r2 != or2 && r2 != NULL) {
            recycleToPool2(tid, r2);
            r2 = NULL;
        }

        if(or1) {
            recycleToPool1(tid, or1);
            or1 = NULL;
        }
        if(or2) {
            recycleToPool2(tid, or2);
            or2 = NULL;
        }
    }

	if(mOptions->split.enabled) {
        // split output by each worker thread
        if(!mOptions->out1.empty())
            config->getWriter1()->writeString(outstr1);
        if(!mOptions->out2.empty())
            config->getWriter2()->writeString(outstr2);
    }

    if(mMergedWriter) {
        // move to heap for writer thread ownership
        mMergedWriter->input(tid, new string(std::move(mergedOutput)));
    }

    if(mFailedWriter) {
        mFailedWriter->input(tid, new string(std::move(failedOut)));
    }

    if(mOverlappedWriter) {
        mOverlappedWriter->input(tid, new string(std::move(overlappedOut)));
    }

    // normal output by left/right writer thread
    if(mRightWriter && mLeftWriter) {
        // write PE - move to heap for writer thread ownership
        mLeftWriter->input(tid, new string(std::move(outstr1)));
        mRightWriter->input(tid, new string(std::move(outstr2)));
    } else if(mLeftWriter) {
        // write singleOutput
        mLeftWriter->input(tid, new string(std::move(singleOutput)));
    }
    // output unpaired reads
    if(mUnpairedLeftWriter && mUnpairedRightWriter) {
        mUnpairedLeftWriter->input(tid, new string(std::move(unpairedOut1)));
        mUnpairedRightWriter->input(tid, new string(std::move(unpairedOut2)));
    } else if(mUnpairedLeftWriter) {
        mUnpairedLeftWriter->input(tid, new string(std::move(unpairedOut1)));
    }

    if(mOptions->split.byFileLines)
        config->markProcessed(readPassed);
    else
        config->markProcessed(leftPack->count);

    if(mOptions->merge.enabled) {
        config->addMergedPairs(mergedCount);
    }

    // stack strings auto-cleanup - no manual delete needed

    delete[] leftPack->data;
    delete[] rightPack->data;
    delete leftPack;
    delete rightPack;

    mPackProcessedCounter++;

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

void PairEndProcessor::readerTask(bool isLeft)
{
    if(mOptions->verbose) {
        if(isLeft)
            loginfo("start to load data of read1");
        else
            loginfo("start to load data of read2");
    }
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    Read** data = new Read*[PACK_SIZE];
    memset(data, 0, sizeof(Read*)*PACK_SIZE);
    FastqReader* reader = NULL;
    if(isLeft) {
        reader = new FastqReader(mOptions->in1, true, mOptions->phred64);
        reader->setReadPool(mLeftReadPool);
    }
    else {
        reader = new FastqReader(mOptions->in2, true, mOptions->phred64);
        reader->setReadPool(mRightReadPool);
    }

    int count=0;
    bool needToBreak = false;
    while(true){
        if(shouldStopReading)
            break;
        Read* read = reader->read();
        if(!read || needToBreak){
            // the last pack
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;

            if(isLeft) {
                mLeftInputLists[mLeftPackReadCounter % mOptions->thread]->produce(pack);
                mLeftPackReadCounter++;
            } else {
                mRightInputLists[mRightPackReadCounter % mOptions->thread]->produce(pack);
                mRightPackReadCounter++;
            }
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
            string msg;
            if(isLeft)
                msg = "Read1: ";
            else
                msg = "Read2: ";
            msg += "loaded " + to_string((lastReported/1000000)) + "M reads";
            loginfo(msg);
        }
        // a full pack
        if(count == PACK_SIZE || needToBreak){
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            
            if(isLeft) {
                mLeftInputLists[mLeftPackReadCounter % mOptions->thread]->produce(pack);
                mLeftPackReadCounter++;
            } else {
                mRightInputLists[mRightPackReadCounter % mOptions->thread]->produce(pack);
                mRightPackReadCounter++;
            }

            //re-initialize data for next pack
            data = new Read*[PACK_SIZE];
            memset(data, 0, sizeof(Read*)*PACK_SIZE);
            // if the processor is far behind this reader, sleep and wait to limit memory usage
            if(isLeft) {
                while(mLeftPackReadCounter - mPackProcessedCounter > PACK_IN_MEM_LIMIT){
                    //cerr<<"sleep"<<endl;
                    slept++;
                    usleep(100);
                }
            } else {
                while(mRightPackReadCounter - mPackProcessedCounter > PACK_IN_MEM_LIMIT){
                    //cerr<<"sleep"<<endl;
                    slept++;
                    usleep(100);
                }
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
                    reader.getBytes(bytesRead, bytesTotal);
                    mOptions->split.size *=  (double)bytesTotal / ((double)bytesRead * (double) mOptions->split.number);
                    if(mOptions->split.size <= 0)
                        mOptions->split.size = 1;
                }
            }*/
        }
    }

    for(int t=0; t<mOptions->thread; t++) {
        if(isLeft)
            mLeftInputLists[t]->setProducerFinished();
        else
            mRightInputLists[t]->setProducerFinished();
    }

    if(mOptions->verbose) {
        if(isLeft) {
            mLeftReaderFinished = true;
            loginfo("Read1: loading completed with " + to_string(mLeftPackReadCounter) + " packs");
        }
        else {
            mRightReaderFinished = true;
            loginfo("Read2: loading completed with " + to_string(mRightPackReadCounter) + " packs");
        }
    }
    

    // if the last data initialized is not used, free it
    if(data != NULL)
        delete[] data;
    if(reader != NULL)
        delete reader;
}

void PairEndProcessor::interleavedReaderTask()
{
    if(mOptions->verbose)
        loginfo("start to load data");
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    Read** dataLeft = new Read*[PACK_SIZE];
    Read** dataRight = new Read*[PACK_SIZE];
    memset(dataLeft, 0, sizeof(Read*)*PACK_SIZE);
    memset(dataRight, 0, sizeof(Read*)*PACK_SIZE);
    FastqReaderPair reader(mOptions->in1, mOptions->in2, true, mOptions->phred64,true);
    int count=0;
    bool needToBreak = false;
    ReadPair* pair = new ReadPair();
    while(true){
        reader.read(pair);
        // TODO: put needToBreak here is just a WAR for resolve some unidentified dead lock issue 
        if(pair->eof() || needToBreak){
            // the last pack
            ReadPack* packLeft = new ReadPack;
            ReadPack* packRight = new ReadPack;
            packLeft->data = dataLeft;
            packRight->data = dataRight;
            packLeft->count = count;
            packRight->count = count;

            mLeftInputLists[mLeftPackReadCounter % mOptions->thread]->produce(packLeft);
            mLeftPackReadCounter++;

            mRightInputLists[mRightPackReadCounter % mOptions->thread]->produce(packRight);
            mRightPackReadCounter++;

            dataLeft = NULL;
            dataRight = NULL;
            break;
        }
        dataLeft[count] = pair->mLeft;
        dataRight[count] = pair->mRight;
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
            ReadPack* packLeft = new ReadPack;
            ReadPack* packRight = new ReadPack;
            packLeft->data = dataLeft;
            packRight->data = dataRight;
            packLeft->count = count;
            packRight->count = count;

            mLeftInputLists[mLeftPackReadCounter % mOptions->thread]->produce(packLeft);
            mLeftPackReadCounter++;

            mRightInputLists[mRightPackReadCounter % mOptions->thread]->produce(packRight);
            mRightPackReadCounter++;

            //re-initialize data for next pack
            dataLeft = new Read*[PACK_SIZE];
            dataRight = new Read*[PACK_SIZE];
            memset(dataLeft, 0, sizeof(Read*)*PACK_SIZE);
            memset(dataRight, 0, sizeof(Read*)*PACK_SIZE);
            // if the consumer is far behind this producer, sleep and wait to limit memory usage
            while(mLeftPackReadCounter - mPackProcessedCounter > PACK_IN_MEM_LIMIT){
                //cerr<<"sleep"<<endl;
                slept++;
                usleep(100);
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

    delete pair;

    for(int t=0; t<mOptions->thread; t++) {
        mLeftInputLists[t]->setProducerFinished();
        mRightInputLists[t]->setProducerFinished();
    }

    if(mOptions->verbose) {
        loginfo("interleaved: loading completed with " + to_string(mLeftPackReadCounter) + " packs");
    }

    mLeftReaderFinished = true;
    mRightReaderFinished = true;

    // if the last data initialized is not used, free it
    if(dataLeft != NULL)
        delete[] dataLeft;
    if(dataRight != NULL)
        delete[] dataRight;
}

void PairEndProcessor::processorTask(ThreadConfig* config)
{
    SingleProducerSingleConsumerList<ReadPack*>* inputLeft = config->getLeftInput();
    SingleProducerSingleConsumerList<ReadPack*>* inputRight = config->getRightInput();
    while(true) {
        if(config->canBeStopped()){
            break;
        }
        while(inputLeft->canBeConsumed() && inputRight->canBeConsumed()) {
            ReadPack* dataLeft = inputLeft->consume();
            ReadPack* dataRight = inputRight->consume();
            processPairEnd(dataLeft, dataRight, config);
        }
        if(inputLeft->isProducerFinished() && !inputLeft->canBeConsumed()) {
            break;
        } else if(inputRight->isProducerFinished() && !inputRight->canBeConsumed()) {
            break;
        } else {
            usleep(100);
        }
    }
    inputLeft->setConsumerFinished();
    inputRight->setConsumerFinished();

    mFinishedThreads++;
    if(mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
        loginfo(msg);
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
        if(mOverlappedWriter)
            mOverlappedWriter->setInputCompleted();
    }
    
    if(mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
        loginfo(msg);
    }
}

void PairEndProcessor::writerTask(WriterThread* config)
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
