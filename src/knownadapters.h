#ifndef KNOWN_ADAPTERS_H
#define KNOWN_ADAPTERS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>

inline map<string, string> getKnownAdapter() {
	map<string, string> knownAdapters;

	knownAdapters["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"] = "Illumina TruSeq Adapter Read 1";
	knownAdapters["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"] = "Illumina TruSeq Adapter Read 2";
	
	return knownAdapters;
}
#endif