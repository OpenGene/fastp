#ifndef FASTP_SIMD_H
#define FASTP_SIMD_H

#include <cstddef>

namespace fastp_simd {

// Count quality metrics for a read in one pass.
// qualstr/seqstr: quality and sequence strings of length len.
// qualThreshold: phred+33 encoded quality threshold for "low quality".
// Outputs: lowQualNum, nBaseNum, totalQual (sum of qual-33 values).
void countQualityMetrics(const char* qualstr, const char* seqstr, int len,
                         char qualThreshold, int& lowQualNum, int& nBaseNum,
                         int& totalQual);

// Reverse complement a DNA sequence.
// src: input sequence of length len.
// dst: output buffer of at least len bytes (may alias src).
void reverseComplement(const char* src, char* dst, int len);

// Count adjacent-base differences for low complexity filter.
// Returns the number of positions where data[i] != data[i+1].
int countAdjacentDiffs(const char* data, int len);

// Count mismatches between two byte strings.
// Returns the number of positions where a[i] != b[i], up to len bytes.
int countMismatches(const char* a, const char* b, int len);

}  // namespace fastp_simd

#endif  // FASTP_SIMD_H
