#ifndef PROCESSOR_H
#define PROCESSOR_H

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

class Processor {
public:
  Processor(Options *opt);
  ~Processor();
  bool process();

private:
  Options *mOptions;
};

#endif