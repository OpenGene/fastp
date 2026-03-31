#ifndef FASTQ_CHUNK_INDEX_H
#define FASTQ_CHUNK_INDEX_H

#include <vector>
#include <string>
#include <cstddef>

// Scans an uncompressed FASTQ file and builds a pack offset index.
// Each entry marks the byte offset where a pack of PACK_SIZE reads begins.
// The last pack may contain fewer than PACK_SIZE reads.
class FastqChunkIndex {
public:
    // Scan the file, recording a byte offset every (4 * packSize) newlines.
    // Returns false if the file cannot be opened or is not a regular file.
    bool build(const std::string& filename, int packSize);

    // Number of packs
    size_t packCount() const { return mOffsets.size() > 0 ? mOffsets.size() - 1 : 0; }

    // Byte offset of pack i
    size_t packStart(size_t i) const { return mOffsets[i]; }

    // Byte length of pack i
    size_t packLength(size_t i) const { return mOffsets[i + 1] - mOffsets[i]; }

    // File descriptor (caller must not close it; owned by this object)
    int fd() const { return mFd; }

    // Total file size
    size_t fileSize() const { return mFileSize; }

    // Clean up
    ~FastqChunkIndex();

private:
    std::vector<size_t> mOffsets;  // packCount+1 entries: [0]=0, [last]=fileSize
    int mFd = -1;
    size_t mFileSize = 0;
};

#endif
