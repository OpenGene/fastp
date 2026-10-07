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
        umi = r1->mSeq->substr(0, min(r1->length(), mOptions->umi.length));
        r1->trimFront(umi.length() + mOptions->umi.skip);
    }
    else if(mOptions->umi.location == UMI_LOC_READ2 && r2){
        umi = r2->mSeq->substr(0, min(r2->length(), mOptions->umi.length));
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
        string umi1 = r1->mSeq->substr(0, min(r1->length(), mOptions->umi.length));
        string umiMerged = umi1;
        r1->trimFront(umi1.length() + mOptions->umi.skip);
        if(r2){
            string umi2 = r2->mSeq->substr(0, min(r2->length(), mOptions->umi.length));
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
    if(!mOptions->umi.umiTag.empty()) {
        // SAM-style tag mode, e.g. --umi_tag=XR produces " XR:Z:UMI"
        tag = " " + mOptions->umi.umiTag + ":Z:";
        if(mOptions->umi.prefix.empty())
            tag = tag + umi;
        else
            tag = tag + mOptions->umi.prefix + "_" + umi;
    } else {
        // legacy mode: append the UMI to the read name with the configured delimiter
        string delimiter = mOptions->umi.delimiter;
        if(mOptions->umi.prefix.empty())
            tag = delimiter + umi;
        else
            tag = delimiter + mOptions->umi.prefix + "_" + umi;
    }
    int spacePos = -1;
    for(int i=0; i<r->mName->length(); i++) {
        if(r->mName->at(i) == ' ') {
            spacePos = i;
            break;
        }
    }
    if(spacePos == -1) {
        r->mName->append(tag);
    } else {
        r->mName->insert(spacePos, tag);
    }

}


bool UmiProcessor::test() {
    Options opt;
    UmiProcessor umiProc(&opt);
    bool passed = true;

    // legacy mode: UMI appended to the read name with the default delimiter ':'
    {
        Read* r = new Read(new string("READNAME"), new string("ACGT"), new string("+"), new string("IIII"));
        umiProc.addUmiToName(r, "AATTCCGG");
        passed &= (*r->mName == "READNAME:AATTCCGG");
        delete r;
    }

    // legacy mode: custom delimiter
    {
        opt.umi.delimiter = "|";
        Read* r = new Read(new string("READNAME"), new string("ACGT"), new string("+"), new string("IIII"));
        umiProc.addUmiToName(r, "AATTCCGG");
        passed &= (*r->mName == "READNAME|AATTCCGG");
        delete r;
    }

    // tag mode: UMI written as a SAM-style tag after the read name
    {
        opt.umi.umiTag = "XR";
        opt.umi.delimiter = ":"; // delimiter should not affect tag mode
        Read* r = new Read(new string("READNAME"), new string("ACGT"), new string("+"), new string("IIII"));
        umiProc.addUmiToName(r, "TCGACC_GCGTAA");
        passed &= (*r->mName == "READNAME XR:Z:TCGACC_GCGTAA");
        delete r;
    }

    // tag mode: read name already contains fields, tag is inserted after the first field
    {
        Read* r = new Read(new string("READNAME/1 1:N:0:INDEX"), new string("ACGT"), new string("+"), new string("IIII"));
        umiProc.addUmiToName(r, "TCGACC_GCGTAA");
        passed &= (*r->mName == "READNAME/1 XR:Z:TCGACC_GCGTAA 1:N:0:INDEX");
        delete r;
    }

    // tag mode: with prefix, the prefix is kept inside the tag value
    {
        opt.umi.prefix = "UMI";
        Read* r = new Read(new string("READNAME"), new string("ACGT"), new string("+"), new string("IIII"));
        umiProc.addUmiToName(r, "AATTCCGG");
        passed &= (*r->mName == "READNAME XR:Z:UMI_AATTCCGG");
        delete r;
    }

    // tag mode: single UMI (read1 only), the underscore-merged value is kept as-is
    {
        opt.umi.prefix = "";
        Read* r = new Read(new string("READNAME"), new string("ACGT"), new string("+"), new string("IIII"));
        umiProc.addUmiToName(r, "TCGACC");
        passed &= (*r->mName == "READNAME XR:Z:TCGACC");
        delete r;
    }

    return passed;
}
