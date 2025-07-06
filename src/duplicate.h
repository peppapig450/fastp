#ifndef DUPLICATE_H
#define DUPLICATE_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include "read.h"
#include "options.h"
#include "common.h"
#include <atomic>
#include <memory>
#include <array>
#include "prime_table.h"

using namespace std;

constexpr uint64 BASE_BUF_BYTES = 1ULL << 29; // 512MB per buffer
constexpr uint32 BASE_BUF_NUM   = 2;          // number of buffers by default

class Duplicate{
public:
    Duplicate(Options* opt);
    ~Duplicate();

    bool checkRead(Read* r1);
    bool checkPair(Read* r1, Read* r2);

    double getDupRate();

private:
    void seq2intvector(const char* data, int len, uint64* output, int posOffset = 0);
    bool applyBloomFilter(uint64* positions);

private:
    Options* mOptions;
    uint64 mBufLenInBits;
    uint64 mBufLenInBytes;
    uint32 mBufNum;
    struct FreeDeleter {
        void operator()(atomic_uchar* p) const noexcept { free(p); }
    };
    std::unique_ptr<atomic_uchar, FreeDeleter> mDupBuf;
    std::unique_ptr<uint64[]> mPrimeArrays;
    std::array<uint64, MAX_BUF_NUM> mPositions;
    atomic_ulong mTotalReads;
    atomic_ulong mDupReads;
    uint64 mOffsetMask;
    
};

#endif
