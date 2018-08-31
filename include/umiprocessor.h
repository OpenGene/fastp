#ifndef UMI_PROCESSOR_H
#define UMI_PROCESSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "options.h"
#include "read.h"

using namespace std;

class UmiProcessor{
public:
    UmiProcessor(Options* opt);
    ~UmiProcessor();
    void process(Read* r1, Read* r2 = NULL);
    void addUmiToName(Read* r, string umi);
    static bool test();

private:
    Options* mOptions;
};


#endif