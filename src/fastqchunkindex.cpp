#include "fastqchunkindex.h"
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <cstring>

// Scan buffer size: 4MB
static const size_t SCAN_BUF_SIZE = 1 << 22;

bool FastqChunkIndex::build(const std::string& filename, int packSize) {
    mFd = open(filename.c_str(), O_RDONLY);
    if (mFd < 0)
        return false;

    struct stat st;
    if (fstat(mFd, &st) != 0 || !S_ISREG(st.st_mode)) {
        ::close(mFd);
        mFd = -1;
        return false;
    }
    mFileSize = st.st_size;
    if (mFileSize == 0) {
        ::close(mFd);
        mFd = -1;
        return false;
    }

#ifdef __linux__
    posix_fadvise(mFd, 0, mFileSize, POSIX_FADV_SEQUENTIAL);
#endif

    // Every FASTQ record is 4 lines, so one pack = packSize * 4 newlines
    const int newlinesPerPack = packSize * 4;

    mOffsets.clear();
    mOffsets.push_back(0);  // first pack starts at byte 0

    char* buf = new char[SCAN_BUF_SIZE];
    size_t fileOffset = 0;
    int nlCount = 0;

    while (fileOffset < mFileSize) {
        size_t toRead = SCAN_BUF_SIZE;
        if (fileOffset + toRead > mFileSize)
            toRead = mFileSize - fileOffset;

        ssize_t n = pread(mFd, buf, toRead, fileOffset);
        if (n <= 0)
            break;

        // Scan for newlines using memchr for speed
        const char* p = buf;
        const char* end = buf + n;
        while (p < end) {
            const char* nl = (const char*)memchr(p, '\n', end - p);
            if (!nl)
                break;
            nlCount++;
            if (nlCount == newlinesPerPack) {
                // The next pack starts right after this newline
                size_t offset = fileOffset + (nl - buf) + 1;
                if (offset < mFileSize)
                    mOffsets.push_back(offset);
                nlCount = 0;
            }
            p = nl + 1;
        }

        fileOffset += n;
    }

    // Sentinel: the end of the last pack is fileSize
    if (mOffsets.back() != mFileSize)
        mOffsets.push_back(mFileSize);

    delete[] buf;
    return true;
}

FastqChunkIndex::~FastqChunkIndex() {
    if (mFd >= 0) {
        ::close(mFd);
        mFd = -1;
    }
}
