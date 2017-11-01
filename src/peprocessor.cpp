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
}

PairEndProcessor::~PairEndProcessor() {
}


bool PairEndProcessor::process(){

    initPackRepository();
    std::thread producer(std::bind(&PairEndProcessor::producerTask, this));

    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&PairEndProcessor::consumerTask, this));
    }

    producer.join();
    for(int t=0; t<mOptions->thread; t++){
        threads[t]->join();
    }

    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
    }

    return true;
}

bool PairEndProcessor::processPairEnd(ReadPairPack* pack){
    for(int p=0;p<pack->count;p++){
        ReadPair* pair = pack->data[p];
        Read* r1 = pair->mLeft;
        Read* r2 = pair->mRight;

        delete pair;
    }

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

void PairEndProcessor::consumePack(){
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

    processPairEnd(data);


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

void PairEndProcessor::consumerTask()
{
    while(true) {
        std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
        if(mProduceFinished && mRepo.writePos == mRepo.readPos){
            lock.unlock();
            break;
        }
        if(mProduceFinished){
            consumePack();
            lock.unlock();
        } else {
            lock.unlock();
            consumePack();
        }
    }
}
