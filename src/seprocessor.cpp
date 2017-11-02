#include "seprocessor.h"
#include "fastqreader.h"
#include <iostream>
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"

SingleEndProcessor::SingleEndProcessor(Options* opt){
    mOptions = opt;
    mProduceFinished = false;
    mFilter = new Filter(opt);
    mOutStream = NULL;
    mZipFile = NULL;
}

SingleEndProcessor::~SingleEndProcessor() {
    delete mFilter;
}

void SingleEndProcessor::initOutput() {
    if(mOptions->out1.empty())
        return;
    if (FastqReader::isZipFastq(mOptions->out1)){
        mZipFile = gzopen(mOptions->out1.c_str(), "w");
        gzsetparams(mZipFile, mOptions->compression, Z_DEFAULT_STRATEGY);
    }
    else {
        mOutStream = new ofstream();
        mOutStream->open(mOptions->out1.c_str(), ifstream::out);
    }
}

void SingleEndProcessor::closeOutput() {
    if (mZipFile){
        gzflush(mZipFile, Z_FINISH);
        gzclose(mZipFile);
        mZipFile = NULL;
    }
    if (mOutStream) {
        if (mOutStream->is_open()){
            mOutStream->flush();
            mOutStream->close();
        }
        delete mOutStream;
    }
}

void SingleEndProcessor::initConfig(ThreadConfig* config) {
    if(mOptions->out1.empty())
        return;
    if(mOutStream != NULL) {
        config->initWriter(mOutStream);
    } else if(mZipFile != NULL) {
        config->initWriter(mZipFile);
    }
}

bool SingleEndProcessor::process(){

    initOutput();

    initPackRepository();
    std::thread producer(std::bind(&SingleEndProcessor::producerTask, this));

    //TODO: get the correct cycles
    int cycle = 151;
    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        configs[t] = new ThreadConfig(mOptions, cycle, false);
        initConfig(configs[t]);
    }

    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&SingleEndProcessor::consumerTask, this, configs[t]));
    }

    producer.join();
    for(int t=0; t<mOptions->thread; t++){
        threads[t]->join();
    }

    // merge stats
    vector<Stats*> preStats;
    vector<Stats*> postStats;
    for(int t=0; t<mOptions->thread; t++){
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
    }
    Stats* finalPreStats = Stats::merge(preStats);
    Stats* finalPostStats = Stats::merge(postStats);

    cout << "pre filtering stats:"<<endl;
    finalPreStats->print();
    cout << "post filtering stats:"<<endl;
    finalPostStats->print();

    // clean up
    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }

    delete threads;
    delete configs;

    closeOutput();

    return true;
}

bool SingleEndProcessor::processSingleEnd(ReadPack* pack, ThreadConfig* config){
    string outstr;
    for(int p=0;p<pack->count;p++){
        Read* r1 = pack->data[p];

        int lowQualNum = 0;
        int nBaseNum = 0;

        // stats the read before filtering
        config->getPreStats1()->statRead(r1, lowQualNum, nBaseNum, mOptions->qualfilter.qualifiedQual);
        if( mFilter->passFilter(r1, lowQualNum, nBaseNum) ) {
            outstr += r1->toString();

            // stats the read after filtering
            config->getPostStats1()->statRead(r1, lowQualNum, nBaseNum, mOptions->qualfilter.qualifiedQual);
        }
    }
    mOutputMtx.lock();
    if(!mOptions->out1.empty())
        config->getWriter1()->writeString(outstr);
    mOutputMtx.unlock();

    delete pack->data;
    delete pack;

    return true;
}

void SingleEndProcessor::initPackRepository() {
    mRepo.packBuffer = new ReadPack*[PACK_NUM_LIMIT];
    memset(mRepo.packBuffer, 0, sizeof(ReadPack*)*PACK_NUM_LIMIT);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
    mRepo.readCounter = 0;
    
}

void SingleEndProcessor::destroyPackRepository() {
    delete mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void SingleEndProcessor::producePack(ReadPack* pack){
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

void SingleEndProcessor::consumePack(ThreadConfig* config){
    ReadPack* data;
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

    processSingleEnd(data, config);


    if (mRepo.readPos >= PACK_NUM_LIMIT)
        mRepo.readPos = 0;

    mRepo.repoNotFull.notify_all();
}

void SingleEndProcessor::producerTask()
{
    int slept = 0;
    Read** data = new Read*[PACK_SIZE];
    memset(data, 0, sizeof(Read*)*PACK_SIZE);
    FastqReader reader1(mOptions->in1);
    int count=0;
    while(true){
        Read* read = reader1.read();
        if(!read){
            // the last pack
            ReadPack* pack = new ReadPack;
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
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            //re-initialize data for next pack
            data = new Read*[PACK_SIZE];
            memset(data, 0, sizeof(Read*)*PACK_SIZE);
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

void SingleEndProcessor::consumerTask(ThreadConfig* config)
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