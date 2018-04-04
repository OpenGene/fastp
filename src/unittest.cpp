#include "unittest.h"
#include "sequence.h"
#include "fastqreader.h"
#include "read.h"
#include "overlapanalysis.h"
#include "filter.h"
#include "adaptertrimmer.h"
#include "basecorrector.h"
#include "polyx.h"
#include <time.h>

UnitTest::UnitTest(){

}

/**
 * Runs all tests and returns 1 on failure, 0 otherwise
 */
int UnitTest::run(){
    bool passed = true;
    passed &= report(Sequence::test(), "Sequence::test");
    passed &= report(FastqReader::test(), "FastqReader::test");
    passed &= report(Read::test(), "Read::test");
    passed &= report(OverlapAnalysis::test(), "OverlapAnalysis::test");
    passed &= report(Filter::test(), "Filter::test");
    passed &= report(AdapterTrimmer::test(), "AdapterTrimmer::test");
    passed &= report(BaseCorrector::test(), "BaseCorrector::test");
    passed &= report(PolyX::test(), "PolyX::test");
    printf("\n==========================\n");
    printf("%s\n\n", passed?"ALL PASSED":"FAILED");
    
    return (passed ? 0 : 1);
}

bool UnitTest::report(bool result, string message) {
    printf("%s:%s\n\n", message.c_str(), result?" PASSED":" FAILED");
    return result;
}
