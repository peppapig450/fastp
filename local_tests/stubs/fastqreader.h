// Stub FastqReader used for local tests. Includes a special mode for simulating
// a large FASTQ file so that Evaluator::evaluateReadNum can exercise its
// byte-based read number estimation.

#ifndef FASTQREADER_STUB_H
#define FASTQREADER_STUB_H

#include <cstddef>
#include <fstream>
#include <string>

class Read;
class ReadPool;

// Sentinel filename that triggers generation of synthetic reads instead of
// reading from disk.
constexpr const char* MOCK_LARGE_DATASET = "MOCK_LARGE_DATASET";

class FastqReader {
public:
    FastqReader(std::string filename, bool hasQuality = true, bool phred64 = false);
    ~FastqReader();

    Read* read();
    void getBytes(size_t& bytesRead, size_t& bytesTotal);
    std::size_t getFileSize() const noexcept;
    bool eof();
    bool hasNoLineBreakAtEnd();
    void setReadPool(ReadPool* rp);
    static bool isZipFastq(std::string filename);
    static bool isFastq(std::string filename);
    static bool test();

private:
    std::ifstream mIn;
    size_t mBytesRead;
    size_t mFileSize;

    // Fields used when simulating a large dataset
    bool  mSimulateLarge;
    long  mTotalReads;
    long  mGenerated;
    size_t mReadLen;
    size_t mBytesPerRead;
};

#endif
