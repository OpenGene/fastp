// SIMD-accelerated functions for fastp using Google Highway.
// This file is compiled once per SIMD target via foreach_target.h.

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/simd.cpp"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace fastp_simd {
namespace HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

void CountQualityMetricsImpl(const char* qualstr, const char* seqstr, int len,
                             char qualThreshold, int& lowQualNum, int& nBaseNum,
                             int& totalQual) {
    const hn::ScalableTag<uint8_t> d;
    const int N = hn::Lanes(d);

    const auto vThresh = hn::Set(d, static_cast<uint8_t>(qualThreshold));
    const auto vN = hn::Set(d, static_cast<uint8_t>('N'));
    const auto v33 = hn::Set(d, static_cast<uint8_t>(33));

    int lowQual = 0;
    int nBase = 0;
    int qualSum = 0;
    int i = 0;

    // Process in chunks, accumulating into wider types to avoid overflow.
    // Process in blocks of up to 255 vectors to avoid u16 overflow during
    // quality sum accumulation (255 * 255 = 65025 fits in u16).
    const int blockSize = 255 * N;

    for (int blockStart = 0; blockStart < len; blockStart += blockSize) {
        int blockEnd = blockStart + blockSize;
        if (blockEnd > len) blockEnd = len;

        // Use u16 accumulator for quality sum within this block.
        // SumsOf2 pairwise-adds adjacent u8 lanes into u16, avoiding
        // PromoteUpperTo which is unavailable on HWY_SCALAR.
        const hn::ScalableTag<uint16_t> d16;
        auto vQualSum16 = hn::Zero(d16);

        for (i = blockStart; i + N <= blockEnd; i += N) {
            const auto vQual = hn::LoadU(d, reinterpret_cast<const uint8_t*>(qualstr + i));
            const auto vSeq = hn::LoadU(d, reinterpret_cast<const uint8_t*>(seqstr + i));

            // Count bases with quality < threshold
            const auto maskLowQ = hn::Lt(vQual, vThresh);
            lowQual += hn::CountTrue(d, maskLowQ);

            // Count N bases
            const auto maskN = hn::Eq(vSeq, vN);
            nBase += hn::CountTrue(d, maskN);

            // Subtract 33 and accumulate quality into u16 accumulator.
            // SumsOf2 adds adjacent u8 pairs -> u16, works on all targets.
            const auto vQualAdj = hn::Sub(vQual, v33);
            vQualSum16 = hn::Add(vQualSum16, hn::SumsOf2(vQualAdj));
        }

        // Reduce u16 accumulator to scalar via SumsOf2 -> u32 then ReduceSum
        const hn::ScalableTag<uint32_t> d32;
        qualSum += static_cast<int>(
            hn::ReduceSum(d32, hn::SumsOf2(vQualSum16)));
    }

    // Scalar tail
    for (; i < len; i++) {
        uint8_t qual = static_cast<uint8_t>(qualstr[i]);
        qualSum += qual - 33;
        if (qual < static_cast<uint8_t>(qualThreshold)) lowQual++;
        if (seqstr[i] == 'N') nBase++;
    }

    lowQualNum = lowQual;
    nBaseNum = nBase;
    totalQual = qualSum;
}

void ReverseComplementImpl(const char* src, char* dst, int len) {
    const hn::ScalableTag<uint8_t> d;
    const int N = hn::Lanes(d);

    // Complement via single table lookup on low nibble of ASCII code.
    // DNA base low nibbles: A/a=1, C/c=3, T/t=4, G/g=7, N=14.
    // Uppercase and lowercase share the same low nibble, so one table
    // handles both.  Unmapped nibbles default to 'N'.
    HWY_ALIGN static const uint8_t kComplement[16] = {
        'N', 'T', 'N', 'G', 'A', 'N', 'N', 'C',
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
    };
    const auto table = hn::LoadDup128(d, kComplement);
    const auto vMask = hn::Set(d, static_cast<uint8_t>(0x0F));

    int i = 0;
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, reinterpret_cast<const uint8_t*>(src + i));
        const auto indices = hn::And(v, vMask);
        const auto comp = hn::TableLookupBytes(table, indices);

        // Reverse and store at mirrored position
        const auto revComp = hn::Reverse(d, comp);
        hn::StoreU(revComp, d, reinterpret_cast<uint8_t*>(dst + len - i - N));
    }

    // Scalar tail
    for (; i < len; i++) {
        char base = src[i];
        char comp;
        switch (base) {
            case 'A': case 'a': comp = 'T'; break;
            case 'T': case 't': comp = 'A'; break;
            case 'C': case 'c': comp = 'G'; break;
            case 'G': case 'g': comp = 'C'; break;
            default: comp = 'N'; break;
        }
        dst[len - 1 - i] = comp;
    }
}

int CountAdjacentDiffsImpl(const char* data, int len) {
    if (len <= 1) return 0;

    const hn::ScalableTag<uint8_t> d;
    const int N = hn::Lanes(d);

    int diff = 0;
    int i = 0;

    // Compare data[i..i+N-1] with data[i+1..i+N]
    for (; i + N < len; i += N) {
        const auto v1 = hn::LoadU(d, reinterpret_cast<const uint8_t*>(data + i));
        const auto v2 = hn::LoadU(d, reinterpret_cast<const uint8_t*>(data + i + 1));
        const auto maskNe = hn::Ne(v1, v2);
        diff += hn::CountTrue(d, maskNe);
    }

    // Scalar tail
    for (; i < len - 1; i++) {
        if (data[i] != data[i + 1]) diff++;
    }

    return diff;
}

int CountMismatchesImpl(const char* a, const char* b, int len) {
    const hn::ScalableTag<uint8_t> d;
    const int N = hn::Lanes(d);

    int diff = 0;
    int i = 0;

    for (; i + N <= len; i += N) {
        const auto va = hn::LoadU(d, reinterpret_cast<const uint8_t*>(a + i));
        const auto vb = hn::LoadU(d, reinterpret_cast<const uint8_t*>(b + i));
        const auto maskNe = hn::Ne(va, vb);
        diff += hn::CountTrue(d, maskNe);
    }

    // Scalar tail
    for (; i < len; i++) {
        if (a[i] != b[i]) diff++;
    }

    return diff;
}

// NOLINTNEXTLINE(google-readability-namespace-comments)
}  // namespace HWY_NAMESPACE
}  // namespace fastp_simd
HWY_AFTER_NAMESPACE();

// ---- Dynamic dispatch wrappers (compiled once) ----
#if HWY_ONCE

#include <cstdio>
#include <cstring>
#include <cstdint>

namespace fastp_simd {

HWY_EXPORT(CountQualityMetricsImpl);
HWY_EXPORT(ReverseComplementImpl);
HWY_EXPORT(CountAdjacentDiffsImpl);
HWY_EXPORT(CountMismatchesImpl);

void countQualityMetrics(const char* qualstr, const char* seqstr, int len,
                         char qualThreshold, int& lowQualNum, int& nBaseNum,
                         int& totalQual) {
    HWY_DYNAMIC_DISPATCH(CountQualityMetricsImpl)(qualstr, seqstr, len,
                                                   qualThreshold, lowQualNum,
                                                   nBaseNum, totalQual);
}

void reverseComplement(const char* src, char* dst, int len) {
    HWY_DYNAMIC_DISPATCH(ReverseComplementImpl)(src, dst, len);
}

int countAdjacentDiffs(const char* data, int len) {
    return HWY_DYNAMIC_DISPATCH(CountAdjacentDiffsImpl)(data, len);
}

int countMismatches(const char* a, const char* b, int len) {
    return HWY_DYNAMIC_DISPATCH(CountMismatchesImpl)(a, b, len);
}

// ---- Scalar reference implementations for testing ----

static void scalarCountQualityMetrics(const char* qualstr, const char* seqstr,
                                       int len, char qualThreshold,
                                       int& lowQualNum, int& nBaseNum,
                                       int& totalQual) {
    lowQualNum = 0;
    nBaseNum = 0;
    totalQual = 0;
    for (int i = 0; i < len; i++) {
        uint8_t q = static_cast<uint8_t>(qualstr[i]);
        totalQual += q - 33;
        if (q < static_cast<uint8_t>(qualThreshold)) lowQualNum++;
        if (seqstr[i] == 'N') nBaseNum++;
    }
}

static void scalarReverseComplement(const char* src, char* dst, int len) {
    for (int i = 0; i < len; i++) {
        char c;
        switch (src[i]) {
            case 'A': case 'a': c = 'T'; break;
            case 'T': case 't': c = 'A'; break;
            case 'C': case 'c': c = 'G'; break;
            case 'G': case 'g': c = 'C'; break;
            default: c = 'N'; break;
        }
        dst[len - 1 - i] = c;
    }
}

static int scalarCountAdjacentDiffs(const char* data, int len) {
    int diff = 0;
    for (int i = 0; i < len - 1; i++) {
        if (data[i] != data[i + 1]) diff++;
    }
    return diff;
}

static int scalarCountMismatches(const char* a, const char* b, int len) {
    int diff = 0;
    for (int i = 0; i < len; i++) {
        if (a[i] != b[i]) diff++;
    }
    return diff;
}

bool testSimd() {
    bool pass = true;

    // --- reverseComplement ---
    // Basic case
    {
        const char* in = "AAAATTTTCCCCGGGG";
        int len = 16;
        char out[16];
        reverseComplement(in, out, len);
        char ref[16];
        scalarReverseComplement(in, ref, len);
        if (memcmp(out, ref, len) != 0) {
            fprintf(stderr, "FAIL: reverseComplement basic\n");
            pass = false;
        }
    }
    // Mixed case
    {
        const char* in = "AaTtCcGgN";
        int len = 9;
        char out[9];
        reverseComplement(in, out, len);
        char ref[9];
        scalarReverseComplement(in, ref, len);
        if (memcmp(out, ref, len) != 0) {
            fprintf(stderr, "FAIL: reverseComplement mixed case\n");
            pass = false;
        }
    }
    // Empty
    {
        char out[1] = {'X'};
        reverseComplement("", out, 0);
        // Should not crash; out unchanged
    }
    // Length 1
    {
        char out[1];
        reverseComplement("A", out, 1);
        if (out[0] != 'T') {
            fprintf(stderr, "FAIL: reverseComplement len=1\n");
            pass = false;
        }
    }
    // Long string (> typical SIMD width) to exercise SIMD + tail
    {
        const char* in = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
        int len = 68;
        char out[68], ref[68];
        reverseComplement(in, out, len);
        scalarReverseComplement(in, ref, len);
        if (memcmp(out, ref, len) != 0) {
            fprintf(stderr, "FAIL: reverseComplement long string\n");
            pass = false;
        }
    }

    // --- countQualityMetrics ---
    // Basic
    {
        // qual: "IIIII" (phred 40 each), seq: "ACGTN"
        const char* qual = "IIIII";
        const char* seq = "ACGTN";
        int lowQ = 0, nBase = 0, totalQ = 0;
        countQualityMetrics(qual, seq, 5, '5', lowQ, nBase, totalQ);
        int refLowQ = 0, refNBase = 0, refTotalQ = 0;
        scalarCountQualityMetrics(qual, seq, 5, '5', refLowQ, refNBase, refTotalQ);
        if (lowQ != refLowQ || nBase != refNBase || totalQ != refTotalQ) {
            fprintf(stderr, "FAIL: countQualityMetrics basic (lowQ=%d/%d nBase=%d/%d totalQ=%d/%d)\n",
                    lowQ, refLowQ, nBase, refNBase, totalQ, refTotalQ);
            pass = false;
        }
    }
    // All low quality
    {
        const char* qual = "!!!!!!!!!!";  // phred 0
        const char* seq = "AAAAAAAAAA";
        int lowQ = 0, nBase = 0, totalQ = 0;
        countQualityMetrics(qual, seq, 10, '5', lowQ, nBase, totalQ);
        int refLowQ = 0, refNBase = 0, refTotalQ = 0;
        scalarCountQualityMetrics(qual, seq, 10, '5', refLowQ, refNBase, refTotalQ);
        if (lowQ != refLowQ || nBase != refNBase || totalQ != refTotalQ) {
            fprintf(stderr, "FAIL: countQualityMetrics all-low\n");
            pass = false;
        }
    }
    // Long string with mixed content
    {
        const char qual[] = "IIIIII!!!!!IIIII55555NNNNN!!!!!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        const char seq[]  = "ACGTNNACGTNACGTNACGTNACGTNACGTNACGTNACGTNACGTNACGTNACGTNACGTNACGTNACG";
        int len = 68;
        int lowQ = 0, nBase = 0, totalQ = 0;
        countQualityMetrics(qual, seq, len, '5', lowQ, nBase, totalQ);
        int refLowQ = 0, refNBase = 0, refTotalQ = 0;
        scalarCountQualityMetrics(qual, seq, len, '5', refLowQ, refNBase, refTotalQ);
        if (lowQ != refLowQ || nBase != refNBase || totalQ != refTotalQ) {
            fprintf(stderr, "FAIL: countQualityMetrics long (lowQ=%d/%d nBase=%d/%d totalQ=%d/%d)\n",
                    lowQ, refLowQ, nBase, refNBase, totalQ, refTotalQ);
            pass = false;
        }
    }

    // --- countAdjacentDiffs ---
    // All same
    {
        int d = countAdjacentDiffs("AAAAAAAAAA", 10);
        if (d != 0) {
            fprintf(stderr, "FAIL: countAdjacentDiffs all-same got %d\n", d);
            pass = false;
        }
    }
    // All different
    {
        const char* s = "ACACACACAC";
        int d = countAdjacentDiffs(s, 10);
        int ref = scalarCountAdjacentDiffs(s, 10);
        if (d != ref) {
            fprintf(stderr, "FAIL: countAdjacentDiffs all-diff got %d expected %d\n", d, ref);
            pass = false;
        }
    }
    // Length 0 and 1
    {
        if (countAdjacentDiffs("A", 1) != 0) {
            fprintf(stderr, "FAIL: countAdjacentDiffs len=1\n");
            pass = false;
        }
        if (countAdjacentDiffs("", 0) != 0) {
            fprintf(stderr, "FAIL: countAdjacentDiffs len=0\n");
            pass = false;
        }
    }
    // Long string
    {
        const char* s = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
        int len = 68;
        int d = countAdjacentDiffs(s, len);
        int ref = scalarCountAdjacentDiffs(s, len);
        if (d != ref) {
            fprintf(stderr, "FAIL: countAdjacentDiffs long got %d expected %d\n", d, ref);
            pass = false;
        }
    }

    // --- countMismatches ---
    // Identical
    {
        const char* s = "ACGTACGTACGT";
        int d = countMismatches(s, s, 12);
        if (d != 0) {
            fprintf(stderr, "FAIL: countMismatches identical got %d\n", d);
            pass = false;
        }
    }
    // All different
    {
        int d = countMismatches("AAAA", "TTTT", 4);
        if (d != 4) {
            fprintf(stderr, "FAIL: countMismatches all-diff got %d\n", d);
            pass = false;
        }
    }
    // Long string
    {
        const char* a = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
        const char* b = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
        int d = countMismatches(a, b, 68);
        if (d != 0) {
            fprintf(stderr, "FAIL: countMismatches long-identical got %d\n", d);
            pass = false;
        }
    }
    {
        const char* a = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        const char* b = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        int d = countMismatches(a, b, 66);
        if (d != 66) {
            fprintf(stderr, "FAIL: countMismatches long-alldiff got %d\n", d);
            pass = false;
        }
    }
    // Length 0
    {
        int d = countMismatches("A", "T", 0);
        if (d != 0) {
            fprintf(stderr, "FAIL: countMismatches len=0 got %d\n", d);
            pass = false;
        }
    }

    return pass;
}

}  // namespace fastp_simd

#endif  // HWY_ONCE
