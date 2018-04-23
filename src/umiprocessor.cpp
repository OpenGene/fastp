#include "umiprocessor.h"

UmiProcessor::UmiProcessor(Options* opt){
    mOptions = opt;
}


UmiProcessor::~UmiProcessor(){
}

void UmiProcessor::process(Read* r1, Read* r2) {
    if(!mOptions->umi.enabled)
        return;

    string umi;
    if(mOptions->umi.location == UMI_LOC_INDEX1)
        umi = r1->firstIndex();
    else if(mOptions->umi.location == UMI_LOC_INDEX2 && r2)
        umi = r2->lastIndex();
    else if(mOptions->umi.location == UMI_LOC_READ1){
        umi = r1->mSeq.mStr.substr(0, min(r1->length(), mOptions->umi.length));
        r1->trimFront(umi.length() + mOptions->umi.skip);
    }
    else if(mOptions->umi.location == UMI_LOC_READ2 && r2){
        umi = r2->mSeq.mStr.substr(0, min(r2->length(), mOptions->umi.length));
        r2->trimFront(umi.length() + mOptions->umi.skip);
    }
    else if(mOptions->umi.location == UMI_LOC_PER_INDEX){
        string umiMerged = r1->firstIndex();
        if(r2) {
            umiMerged = umiMerged + "_" + r2->lastIndex();
        }

        addUmiToName(r1, umiMerged);
        if(r2) {
            addUmiToName(r2, umiMerged);
        }
    }
    else if(mOptions->umi.location == UMI_LOC_PER_READ){
        string umi1 = r1->mSeq.mStr.substr(0, min(r1->length(), mOptions->umi.length));
        string umiMerged = umi1;
        r1->trimFront(umi1.length() + mOptions->umi.skip);
        if(r2){
            string umi2 = r2->mSeq.mStr.substr(0, min(r2->length(), mOptions->umi.length));
            umiMerged = umiMerged + "_" + umi2;
            r2->trimFront(umi2.length() + mOptions->umi.skip);
        }

        addUmiToName(r1, umiMerged);
        if(r2){
            addUmiToName(r2, umiMerged);
        }
    }

    if(mOptions->umi.location != UMI_LOC_PER_INDEX && mOptions->umi.location != UMI_LOC_PER_READ) {
        if(r1 && !umi.empty()) 
            addUmiToName(r1, umi);
        if(r2 && !umi.empty())
            addUmiToName(r2, umi);
    }
}

void UmiProcessor::addUmiToName(Read* r, string umi){
    string tag;
    if(mOptions->umi.prefix.empty())
        tag = ":" + umi;
    else
        tag = ":" + mOptions->umi.prefix + "_" + umi;
    int spacePos = -1;
    for(int i=0; i<r->mName.length(); i++) {
        if(r->mName[i] == ' ') {
            spacePos = i;
            break;
        }
    }
    if(spacePos == -1) {
        r->mName = r->mName + tag;
    } else {
        r->mName = r->mName.substr(0, spacePos) + tag + r->mName.substr(spacePos, r->mName.length() - spacePos);
    }

}


bool UmiProcessor::test() {
    return true;
}