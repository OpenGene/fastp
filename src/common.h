#ifndef COMMON_H
#define COMMON_H

#define FASTP_VER "1.3.7"

#define _DEBUG false

#ifndef _WIN32
	typedef long int64;
	typedef unsigned long uint64;
#else
	typedef long long int64;
	typedef unsigned long long uint64;
#endif

typedef int int32;
typedef unsigned int uint32;

typedef short int16;
typedef unsigned short uint16;

typedef char int8;
typedef unsigned char uint8;

const char ATCG_BASES[] = {'A', 'T', 'C', 'G'};

// how many reads one pack has
// ~1000 reads × ~350 bytes/read ≈ 350KB, near FASTQ LZ77 saturation (~256KB),
// so each pack compresses efficiently as a single gzip member.
static const int PACK_SIZE = 1000;

// if one pack is produced, but not consumed, it will be kept in the memory
// this number limit the number of in memory packs
// if the number of in memory packs is full, the producer thread should sleep
// Scaled down from 128 to match PACK_SIZE increase (256→1000) and maintain
// similar peak memory: 32 × 1000 reads ≈ 32K reads in flight ≈ old 128 × 256.
static const int PACK_IN_MEM_LIMIT = 32;

// Each worker thread owns one input list, and a list's newest item only becomes
// consumable once another item is produced behind it (or the producer finishes).
// Reader backpressure must therefore allow at least one in-flight pack per worker;
// otherwise, with more workers than PACK_IN_MEM_LIMIT, readers stop before any
// worker has a consumable pack and all threads wait forever (#721).
inline long packInMemLimit(int threads) {
    return threads * 2L > PACK_IN_MEM_LIMIT ? threads * 2L : PACK_IN_MEM_LIMIT;
}

// Total buffered reads is roughly packInMemLimit(threads) * packSize(threads).
// packInMemLimit can't safely drop below 2*threads (see above), so it grows
// linearly with thread count with no slack to spare. To keep total buffered
// memory from also growing unboundedly on high-thread-count hosts, shrink the
// per-pack read count instead once thread count exceeds the point where
// packInMemLimit's own floor (PACK_IN_MEM_LIMIT) stops mattering. This keeps
// packInMemLimit(threads) * packSize(threads) roughly constant rather than
// linear in threads, at the cost of smaller (so slightly less efficient) gzip
// members per pack at very high thread counts. PACK_SIZE_FLOOR keeps packs
// from shrinking so far that per-pack overhead dominates.
static const int PACK_SIZE_FLOOR = 100;
inline int packSize(int threads) {
    const int baselineThreads = PACK_IN_MEM_LIMIT / 2;
    if (threads <= baselineThreads)
        return PACK_SIZE;
    long scaled = (long)PACK_SIZE * baselineThreads / threads;
    return scaled > PACK_SIZE_FLOOR ? (int)scaled : PACK_SIZE_FLOOR;
}

// different filtering results, bigger number means worse
// if r1 and r2 are both failed, then the bigger one of the two results will be recorded
// we reserve some gaps for future types to be added
static const int PASS_FILTER = 0;
static const int FAIL_POLY_X = 4;
static const int FAIL_OVERLAP = 8;
static const int FAIL_N_BASE = 12;
static const int FAIL_LENGTH = 16;
static const int FAIL_TOO_LONG = 17;
static const int FAIL_QUALITY = 20;
static const int FAIL_COMPLEXITY = 24;
static const int FAIL_ADAPTER_DIMER = 28;

// how many types in total we support
static const int FILTER_RESULT_TYPES = 32;

const static char* FAILED_TYPES[FILTER_RESULT_TYPES] = {
	"passed", "", "", "",
	"failed_polyx_filter", "", "", "",
	"failed_bad_overlap", "", "", "",
	"failed_too_many_n_bases", "", "", "",
	"failed_too_short", "failed_too_long", "", "",
	"failed_quality_filter", "", "", "",
	"failed_low_complexity", "", "", "",
	"failed_adapter_dimer", "", "", ""
};

#endif /* COMMON_H */
