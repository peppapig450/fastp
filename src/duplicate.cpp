#include "duplicate.h"
#include "overlapanalysis.h"
#include <cstring>
#include <cstdlib>
#include <math.h>
#include <algorithm>
#include <cassert>
#include "util.h"

static constexpr uint64_t A_HASH = 7;
static constexpr uint64_t T_HASH = 222;
static constexpr uint64_t C_HASH = 74;
static constexpr uint64_t G_HASH = 31;
static constexpr uint64_t N_HASH = 13;

inline uint64_t baseHash(char b) {
    switch(b) {
        case 'A': return A_HASH;
        case 'T': return T_HASH;
        case 'C': return C_HASH;
        case 'G': return G_HASH;
        default:  return N_HASH;
    }
}

#include "prime_table.h"

Duplicate::Duplicate(Options* opt) {
    mOptions = opt;

    // 1G mem required
    mBufLenInBytes = BASE_BUF_BYTES;
    mBufNum = BASE_BUF_NUM;

    // the memory usage increases as the accuracy level increases
    // level 1: 1G
    // level 2: 2G
    // level 3: 4G
    // level 4: 8G
    // level 5: 16G
    // level 6: 24G
    switch(mOptions->duplicate.accuracyLevel) {
        case 1:
            break;
        case 2:
            mBufLenInBytes *= 2;
            break;
        case 3:
            mBufLenInBytes *= 2;
            mBufNum *= 2;
            break;
        case 4:
            mBufLenInBytes *= 4;
            mBufNum *= 2;
            break;
        case 5:
            mBufLenInBytes *= 8;
            mBufNum *= 2;
            break;
        case 6:
            mBufLenInBytes *= 8;
            mBufNum *= 3;
            break;
        default:
            break;
    }

    mOffsetMask = PRIME_ARRAY_LEN * mBufNum - 1;

    mBufLenInBits = mBufLenInBytes << 3;
    void* raw = nullptr;
    if(posix_memalign(&raw, 64, mBufLenInBytes * mBufNum)) {
        error_exit("Out of memory, failed to allocate " + to_string(mBufLenInBytes * mBufNum) + " bytes buffer for duplication analysis, please reduce the dup_accuracy_level and try again.");
    }
    mDupBuf.reset(static_cast<atomic_uchar*>(raw));
    for(uint64 i = 0; i < mBufLenInBytes * mBufNum; ++i) {
        mDupBuf.get()[i] = 0;
    }

    mPrimeArrays = std::unique_ptr<uint64[]>(new uint64[mBufNum * PRIME_ARRAY_LEN]);
    std::copy_n(PRIME_TABLE.begin(), mBufNum * PRIME_ARRAY_LEN, mPrimeArrays.get());

    mTotalReads = 0;
    mDupReads = 0;
}


Duplicate::~Duplicate(){
    // resources are managed by smart pointers
}

void Duplicate::seq2intvector(const char* data, int len, uint64* output, int posOffset) {
    for(int p=0; p<len; p++) {
        uint64 base = baseHash(data[p]);
        for(int i=0; i<mBufNum; i++) {
            assert(i < mBufNum);
            int offset = (p+posOffset)*mBufNum + i;
            offset &= mOffsetMask;
            assert(offset < mBufNum * PRIME_ARRAY_LEN);
            output[i] += mPrimeArrays[offset] * (base + (p+posOffset));
        }
    }
}

bool Duplicate::checkRead(Read* r) {
    fill(mPositions.begin(), mPositions.begin() + mBufNum, 0);
    int len = r->length();
    seq2intvector(r->mSeq->c_str(), len, mPositions.data());
    bool isDup = applyBloomFilter(mPositions.data());

    mTotalReads++;
    if(isDup)
        mDupReads++;

    return isDup;
}

bool Duplicate::checkPair(Read* r1, Read* r2) {
    fill(mPositions.begin(), mPositions.begin() + mBufNum, 0);
    seq2intvector(r1->mSeq->c_str(), r1->length(), mPositions.data());
    seq2intvector(r2->mSeq->c_str(), r2->length(), mPositions.data(), r1->length());
    bool isDup = applyBloomFilter(mPositions.data());

    mTotalReads++;
    if(isDup)
        mDupReads++;

    return isDup;
}

bool Duplicate::applyBloomFilter(uint64* positions) {
    bool isDup = true;
    for(int i=0; i<mBufNum; i++) {
        uint64 pos = positions[i] % mBufLenInBits;
        uint64 bytePos = pos >> 3;
        uint32 bitOffset = pos & 0x07;
        uint8 byte = (0x01) << bitOffset;

        //isDup = isDup && (mDupBuf[i * mBufLenInBytes + bytePos] & byte);
        uint8 ret = atomic_fetch_or(mDupBuf.get() + i * mBufLenInBytes + bytePos, byte);
        isDup = (ret & byte) != 0;
    }
    return isDup;
}

double Duplicate::getDupRate() {
    if(mTotalReads == 0)
        return 0.0;
    return (double)mDupReads/(double)mTotalReads;
}
