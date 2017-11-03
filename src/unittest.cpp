#include "unittest.h"
#include "sequence.h"
#include "fastqreader.h"
#include "read.h"
#include "overlapanalysis.h"
#include <time.h>

UnitTest::UnitTest(){

}

void UnitTest::run(){
    bool passed = true;
    passed &= report(Sequence::test(), "Sequence::test");
    passed &= report(FastqReader::test(), "FastqReader::test");
    passed &= report(Read::test(), "Read::test");
    passed &= report(OverlapAnalysis::test(), "OverlapAnalysis::test");
    printf("\n==========================\n");
    printf("%s\n\n", passed?"ALL PASSED":"FAILED");
}

bool UnitTest::report(bool result, string message) {
    printf("%s:%s\n\n", message.c_str(), result?" PASSED":" FAILED");
    return result;
}