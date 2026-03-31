#include "fastqchunkparser.h"
#include "common.h"
#include <cstring>
#include <vector>
#include <iostream>

using namespace std;

// Find next newline in buf[pos..len), return position of char after \n.
// Handles \r\n. Returns len if no newline found.
static inline size_t nextLine(const char* buf, size_t pos, size_t len) {
    const char* p = (const char*)memchr(buf + pos, '\n', len - pos);
    if (!p)
        return len;
    return (p - buf) + 1;
}

// Extract a line from buf[start..end) stripping trailing \r\n
static inline void extractLine(const char* buf, size_t start, size_t lineEnd, string* out) {
    size_t end = lineEnd;
    // strip trailing \n and \r
    while (end > start && (buf[end - 1] == '\n' || buf[end - 1] == '\r'))
        end--;
    out->assign(buf + start, end - start);
}

ReadPack* FastqChunkParser::parse(const char* buf, size_t len, int tid, ReadPool* pool, bool phred64) {
    // Pre-count records: count newlines / 4
    int estimatedReads = 0;
    {
        const char* p = buf;
        const char* end = buf + len;
        while (p < end) {
            p = (const char*)memchr(p, '\n', end - p);
            if (!p) break;
            estimatedReads++;
            p++;
        }
        estimatedReads /= 4;
    }
    if (estimatedReads <= 0)
        estimatedReads = 1;

    Read** data = new Read*[estimatedReads];
    int count = 0;
    size_t pos = 0;

    while (pos < len) {
        // Skip empty lines
        while (pos < len && (buf[pos] == '\n' || buf[pos] == '\r'))
            pos++;
        if (pos >= len)
            break;

        // Line 1: name (starts with @)
        size_t nameStart = pos;
        size_t nameEnd = nextLine(buf, pos, len);
        if (nameStart >= len || buf[nameStart] != '@')
            break;
        pos = nameEnd;

        // Line 2: sequence
        if (pos >= len) break;
        size_t seqStart = pos;
        size_t seqEnd = nextLine(buf, pos, len);
        pos = seqEnd;

        // Line 3: strand (+)
        if (pos >= len) break;
        size_t strandStart = pos;
        size_t strandEnd = nextLine(buf, pos, len);
        pos = strandEnd;

        // Line 4: quality
        if (pos >= len && strandEnd >= len) break;
        size_t qualStart = pos;
        size_t qualEnd = nextLine(buf, pos, len);
        pos = qualEnd;

        // Build Read object
        Read* r = nullptr;
        if (pool)
            r = pool->getOne();

        if (r) {
            extractLine(buf, nameStart, nameEnd, r->mName);
            extractLine(buf, seqStart, seqEnd, r->mSeq);
            extractLine(buf, strandStart, strandEnd, r->mStrand);
            extractLine(buf, qualStart, qualEnd, r->mQuality);
        } else {
            string* name = new string();
            string* seq = new string();
            string* strand = new string();
            string* quality = new string();
            extractLine(buf, nameStart, nameEnd, name);
            extractLine(buf, seqStart, seqEnd, seq);
            extractLine(buf, strandStart, strandEnd, strand);
            extractLine(buf, qualStart, qualEnd, quality);
            r = new Read(name, seq, strand, quality, phred64);
        }

        if (count >= estimatedReads) {
            // Shouldn't happen, but be safe
            Read** newData = new Read*[estimatedReads * 2];
            memcpy(newData, data, sizeof(Read*) * count);
            delete[] data;
            data = newData;
            estimatedReads *= 2;
        }
        data[count++] = r;
    }

    ReadPack* pack = new ReadPack;
    pack->data = data;
    pack->count = count;
    return pack;
}
