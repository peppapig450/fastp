// Stubbed FastqReader implementation with support for generating synthetic
// reads for a virtual large dataset.

#include "fastqreader.h"
#include "../../src/read.h"
#include <sys/stat.h>
#include <string>

FastqReader::FastqReader(std::string filename, bool hasQuality, bool phred64)
    : mIn(), mBytesRead(0), mFileSize(0), mSimulateLarge(false),
      mTotalReads(0), mGenerated(0), mReadLen(0), mBytesPerRead(0) {
    if (filename == MOCK_LARGE_DATASET) {
        mSimulateLarge = true;
        mReadLen = 100;                    // constant read length
        std::string name = "@SEQ";        // constant read name
        std::string seq(mReadLen, 'A');
        std::string plus = "+";
        std::string qual(mReadLen, 'I');
        mBytesPerRead = name.size() + 1 + seq.size() + 1 + plus.size() + 1 +
                        qual.size() + 1;
        mTotalReads = 1024 * 1024;        // 1M reads, exceeds READ_LIMIT
        mFileSize = mTotalReads * mBytesPerRead;
    } else {
        mIn.open(filename);
        mIn.seekg(0, std::ios::end);
        mFileSize = static_cast<size_t>(mIn.tellg());
        mIn.seekg(0);
    }
}

FastqReader::~FastqReader() {}

Read* FastqReader::read() {
    if (mSimulateLarge) {
        if (mGenerated >= mTotalReads) return nullptr;
        mGenerated++;
        mBytesRead += mBytesPerRead;
        std::string name = "@SEQ";
        std::string seq(mReadLen, 'A');
        std::string qual(mReadLen, 'I');
        return new Read(name, seq, qual);
    }

    if (!mIn.good()) return nullptr;
    std::string name, seq, plus, qual;
    if (!std::getline(mIn, name)) return nullptr;
    if (!std::getline(mIn, seq)) return nullptr;
    if (!std::getline(mIn, plus)) return nullptr;
    if (!std::getline(mIn, qual)) return nullptr;
    mBytesRead = static_cast<size_t>(mIn.tellg());
    return new Read(name, seq, qual);
}

void FastqReader::getBytes(size_t& bytesRead, size_t& bytesTotal) {
    bytesRead = mBytesRead;
    bytesTotal = mFileSize;
}

std::size_t FastqReader::getFileSize() const noexcept {
    return mFileSize;
}

bool FastqReader::eof() {
    if (mSimulateLarge) return mGenerated >= mTotalReads;
    return mIn.eof();
}

bool FastqReader::hasNoLineBreakAtEnd() { return false; }

void FastqReader::setReadPool(ReadPool* rp) {}

bool FastqReader::isZipFastq(std::string filename) { return false; }

bool FastqReader::isFastq(std::string filename) { return true; }

bool FastqReader::test() { return true; }

