#ifndef JSON_REPORTER_H
#define JSON_REPORTER_H

#include "common.h"
#include "filterresult.h"
#include "options.h"
#include "stats.h"
#include "util.h"
#include <atomic>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

class JsonReporter {
public:
  JsonReporter(Options *opt);
  ~JsonReporter();

  void setDup(double dupRate);
  void setInsertHist(atomic_long *insertHist, int insertSizePeak);
  void report(FilterResult *result, Stats *preStats1, Stats *postStats1,
              Stats *preStats2 = NULL, Stats *postStats2 = NULL);

private:
  Options *mOptions;
  int *mDupHist;
  double *mDupMeanGC;
  double mDupRate;
  atomic_long *mInsertHist;
  int mInsertSizePeak;
};

#endif
