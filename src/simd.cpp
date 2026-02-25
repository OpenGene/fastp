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

        // Use u16 accumulators for quality sum within this block
        const hn::ScalableTag<uint16_t> d16;
        auto vQualSumLo = hn::Zero(d16);
        auto vQualSumHi = hn::Zero(d16);

        for (i = blockStart; i + N <= blockEnd; i += N) {
            const auto vQual = hn::Load(d, reinterpret_cast<const uint8_t*>(qualstr + i));
            const auto vSeq = hn::Load(d, reinterpret_cast<const uint8_t*>(seqstr + i));

            // Count bases with quality < threshold
            const auto maskLowQ = hn::Lt(vQual, vThresh);
            lowQual += hn::CountTrue(d, maskLowQ);

            // Count N bases
            const auto maskN = hn::Eq(vSeq, vN);
            nBase += hn::CountTrue(d, maskN);

            // Subtract 33 and accumulate quality into u16 accumulators
            const auto vQualAdj = hn::Sub(vQual, v33);
            // Widen to u16 and add: split into lower and upper halves
            const auto lo = hn::PromoteLowerTo(d16, vQualAdj);
            const auto hi = hn::PromoteUpperTo(d16, vQualAdj);
            vQualSumLo = hn::Add(vQualSumLo, lo);
            vQualSumHi = hn::Add(vQualSumHi, hi);
        }

        // Reduce u16 accumulators to scalar
        // Widen to u32 first to avoid overflow in ReduceSum
        const hn::ScalableTag<uint32_t> d32;
        qualSum += static_cast<int>(
            hn::ReduceSum(d32, hn::PromoteLowerTo(d32, vQualSumLo)) +
            hn::ReduceSum(d32, hn::PromoteUpperTo(d32, vQualSumLo)) +
            hn::ReduceSum(d32, hn::PromoteLowerTo(d32, vQualSumHi)) +
            hn::ReduceSum(d32, hn::PromoteUpperTo(d32, vQualSumHi)));
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

    // Complement via comparisons — portable across all SIMD widths
    const auto vA = hn::Set(d, static_cast<uint8_t>('A'));
    const auto vT = hn::Set(d, static_cast<uint8_t>('T'));
    const auto vC = hn::Set(d, static_cast<uint8_t>('C'));
    const auto vG = hn::Set(d, static_cast<uint8_t>('G'));
    const auto va = hn::Set(d, static_cast<uint8_t>('a'));
    const auto vt = hn::Set(d, static_cast<uint8_t>('t'));
    const auto vc = hn::Set(d, static_cast<uint8_t>('c'));
    const auto vg = hn::Set(d, static_cast<uint8_t>('g'));
    const auto vN = hn::Set(d, static_cast<uint8_t>('N'));

    int i = 0;
    for (; i + N <= len; i += N) {
        const auto v = hn::Load(d, reinterpret_cast<const uint8_t*>(src + i));

        // Build complement using conditional selects
        auto comp = vN;  // default: N
        comp = hn::IfThenElse(hn::Eq(v, vA), vT, comp);
        comp = hn::IfThenElse(hn::Eq(v, vT), vA, comp);
        comp = hn::IfThenElse(hn::Eq(v, vC), vG, comp);
        comp = hn::IfThenElse(hn::Eq(v, vG), vC, comp);
        comp = hn::IfThenElse(hn::Eq(v, va), vT, comp);
        comp = hn::IfThenElse(hn::Eq(v, vt), vA, comp);
        comp = hn::IfThenElse(hn::Eq(v, vc), vG, comp);
        comp = hn::IfThenElse(hn::Eq(v, vg), vC, comp);

        // Reverse and store at mirrored position
        const auto revComp = hn::Reverse(d, comp);
        hn::Store(revComp, d, reinterpret_cast<uint8_t*>(dst + len - i - N));
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
        const auto v1 = hn::Load(d, reinterpret_cast<const uint8_t*>(data + i));
        const auto v2 = hn::Load(d, reinterpret_cast<const uint8_t*>(data + i + 1));
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
        const auto va = hn::Load(d, reinterpret_cast<const uint8_t*>(a + i));
        const auto vb = hn::Load(d, reinterpret_cast<const uint8_t*>(b + i));
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

}  // namespace fastp_simd

#endif  // HWY_ONCE
