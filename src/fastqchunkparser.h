#ifndef FASTQ_CHUNK_PARSER_H
#define FASTQ_CHUNK_PARSER_H

#include "read.h"
#include "readpool.h"

// Parse a raw byte buffer (from pread) into a ReadPack.
// The buffer must start at a record boundary and contain complete records.
class FastqChunkParser {
public:
    // Parse buf[0..len) into Read* objects, return a ReadPack.
    // If pool is non-null, try to reuse Read objects from the pool.
    // tid is the thread id for the pool.
    static ReadPack* parse(const char* buf, size_t len, int tid, ReadPool* pool, bool phred64);
};

#endif
