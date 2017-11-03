#include "peprocessor.h"
#include "fastqreader.h"
#include <iostream>
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"

PairEndProcessor::PairEndProcessor(Options* opt){
    mOptions = opt;
    mProduceFinished = false;
    mFilter = new Filter(opt);
    mOutStream1 = NULL;
    mZipFile1 = NULL;
    mOutStream2 = NULL;
    mZipFile2 = NULL;
}

PairEndProcessor::~PairEndProcessor() {
}

void PairEndProcessor::initOutput() {
    if(mOptions->out1.empty() || mOptions->out2.empty())
        return;
    if (FastqReader::isZipFastq(mOptions->out1)){
        mZipFile1 = gzopen(mOptions->out1.c_str(), "w");
        gzsetparams(mZipFile1, mOptions->compression, Z_DEFAULT_STRATEGY);
        gzbuffer(mZipFile1, 1024*1024);
        mZipFile2 = gzopen(mOptions->out2.c_str(), "w");
        gzsetparams(mZipFile2, mOptions->compression, Z_DEFAULT_STRATEGY);
        gzbuffer(mZipFile2, 1024*1024);
    }
    else {
        mOutStream1 = new ofstream();
        mOutStream1->open(mOptions->out1.c_str(), ifstream::out);
        mOutStream2 = new ofstream();
        mOutStream2->open(mOptions->out2.c_str(), ifstream::out);
    }
}

void PairEndProcessor::closeOutput() {
    if (mZipFile1){
        gzflush(mZipFile1, Z_FINISH);
        gzclose(mZipFile1);
        mZipFile1 = NULL;
    }
    if (mZipFile2){
        gzflush(mZipFile2, Z_FINISH);
        gzclose(mZipFile2);
        mZipFile2 = NULL;
    }
    if (mOutStream1) {
        if (mOutStream1->is_open()){
            mOutStream1->flush();
            mOutStream1->close();
        }
        delete mOutStream1;
    }
    if (mOutStream2) {
        if (mOutStream2->is_open()){
            mOutStream2->flush();
            mOutStream2->close();
        }
        delete mOutStream2;
    }
}

void PairEndProcessor::initConfig(ThreadConfig* config) {
    if(mOptions->out1.empty())
        return;
    if(mOutStream1 != NULL && mOutStream2 != NULL) {
        config->initWriter(mOutStream1, mOutStream2);
    } else if(mZipFile1 != NULL && mZipFile2 != NULL) {
        config->initWriter(mZipFile1, mZipFile2);
    }
}


bool PairEndProcessor::process(){
    initOutput();

    initPackRepository();
    std::thread producer(std::bind(&PairEndProcessor::producerTask, this));

    //TODO: get the correct cycles
    int cycle = 151;
    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        configs[t] = new ThreadConfig(mOptions, cycle, true);
        initConfig(configs[t]);
    }

    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&PairEndProcessor::consumerTask, this, configs[t]));
    }

    producer.join();
    for(int t=0; t<mOptions->thread; t++){
        threads[t]->join();
    }

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

    cout << "pre filtering stats1:"<<endl;
    finalPreStats1->print();
    cout << endl;
    cout << "post filtering stats1:"<<endl;
    finalPostStats1->print();
    cout << endl;
    cout << "pre filtering stats2:"<<endl;
    finalPreStats2->print();
    cout << endl;
    cout << "post filtering stats2:"<<endl;
    finalPostStats2->print();

    cout << endl;
    cout << "filter results:"<<endl;
    finalFilterResult->print();

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

    delete threads;
    delete configs;

    closeOutput();

    return true;
}

bool PairEndProcessor::processPairEnd(ReadPairPack* pack, ThreadConfig* config){
    string outstr1;
    string outstr2;
    for(int p=0;p<pack->count;p++){
        ReadPair* pair = pack->data[p];
        Read* or1 = pair->mLeft;
        Read* or2 = pair->mRight;

        int lowQualNum1 = 0;
        int nBaseNum1 = 0;
        int lowQualNum2 = 0;
        int nBaseNum2 = 0;

        // stats the original read before trimming
        config->getPreStats1()->statRead(or1, lowQualNum1, nBaseNum1, mOptions->qualfilter.qualifiedQual);
        config->getPreStats2()->statRead(or2, lowQualNum2, nBaseNum2, mOptions->qualfilter.qualifiedQual);

        // trim in head and tail, and cut adapters
        Read* r1 = mFilter->trimAndCutAdapter(or1);
        Read* r2 = mFilter->trimAndCutAdapter(or2);

        int result1 = mFilter->passFilter(r1, lowQualNum1, nBaseNum1);
        int result2 = mFilter->passFilter(r2, lowQualNum2, nBaseNum2);

        config->addFilterResult(max(result1, result2));

        if( r1 != NULL &&  result1 == PASS_FILTER && r2 != NULL && result2 == PASS_FILTER ) {
            
            outstr1 += r1->toString();
            outstr2 += r2->toString();

            // stats the read after filtering
            config->getPostStats1()->statRead(r1, lowQualNum1, nBaseNum1, mOptions->qualfilter.qualifiedQual);
            config->getPostStats2()->statRead(r2, lowQualNum2, nBaseNum2, mOptions->qualfilter.qualifiedQual);
        }


        delete pair;
        // if no trimming applied, r1 should be identical to or1
        if(r1 != or1 && r1 != NULL)
            delete r1;
        // if no trimming applied, r1 should be identical to or1
        if(r2 != or2 && r2 != NULL)
            delete r2;
    }

    mOutputMtx.lock();
    if(!mOptions->out1.empty())
        config->getWriter1()->writeString(outstr1);
    if(!mOptions->out2.empty())
        config->getWriter2()->writeString(outstr2);
    mOutputMtx.unlock();

    delete pack->data;
    delete pack;

    return true;
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
    mRepo.readCounter = 0;
    
}

void PairEndProcessor::destroyPackRepository() {
    delete mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void PairEndProcessor::producePack(ReadPairPack* pack){
    std::unique_lock<std::mutex> lock(mRepo.mtx);
    while(((mRepo.writePos + 1) % PACK_NUM_LIMIT)
        == mRepo.readPos) {
        mRepo.repoNotFull.wait(lock);
    }

    mRepo.packBuffer[mRepo.writePos] = pack;
    mRepo.writePos++;

    if (mRepo.writePos == PACK_NUM_LIMIT)
        mRepo.writePos = 0;

    mRepo.repoNotEmpty.notify_all();
    lock.unlock();
}

void PairEndProcessor::consumePack(ThreadConfig* config){
    ReadPairPack* data;
    std::unique_lock<std::mutex> lock(mRepo.mtx);
    // read buffer is empty, just wait here.
    while(mRepo.writePos == mRepo.readPos) {
        if(mProduceFinished){
            lock.unlock();
            return;
        }
        mRepo.repoNotEmpty.wait(lock);
    }

    data = mRepo.packBuffer[mRepo.readPos];
    (mRepo.readPos)++;
    lock.unlock();

    processPairEnd(data, config);


    if (mRepo.readPos >= PACK_NUM_LIMIT)
        mRepo.readPos = 0;

    mRepo.repoNotFull.notify_all();
}

void PairEndProcessor::producerTask()
{
    int slept = 0;
    ReadPair** data = new ReadPair*[PACK_SIZE];
    memset(data, 0, sizeof(ReadPair*)*PACK_SIZE);
    FastqReaderPair reader(mOptions->in1, mOptions->in2);
    int count=0;
    while(true){
        ReadPair* read = reader.read();
        if(!read){
            // the last pack
            ReadPairPack* pack = new ReadPairPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            data = NULL;
            break;
        }
        data[count] = read;
        count++;
        // a full pack
        if(count == PACK_SIZE){
            ReadPairPack* pack = new ReadPairPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            //re-initialize data for next pack
            data = new ReadPair*[PACK_SIZE];
            memset(data, 0, sizeof(ReadPair*)*PACK_SIZE);
            // reset count to 0
            count = 0;
            // if the consumer is far behind this producer, sleep and wait to limit memory usage
            while(mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT){
                //cout<<"sleep"<<endl;
                slept++;
                usleep(1000);
            }
        }
    }

    std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
    mProduceFinished = true;
    lock.unlock();

    // if the last data initialized is not used, free it
    if(data != NULL)
        delete data;
}

void PairEndProcessor::consumerTask(ThreadConfig* config)
{
    while(true) {
        std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
        if(mProduceFinished && mRepo.writePos == mRepo.readPos){
            lock.unlock();
            break;
        }
        if(mProduceFinished){
            consumePack(config);
            lock.unlock();
        } else {
            lock.unlock();
            consumePack(config);
        }
    }
}
