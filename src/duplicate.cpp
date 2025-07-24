#include "duplicate.h"
#include "overlapanalysis.h"
#include <memory.h>
#include <math.h>
#include "util.h"

const int PRIME_ARRAY_LEN = 1<<9;

#ifdef DEBUG_DUPLICATE
void runAllDebugTests(Duplicate* dup);  // Forward declaration
#endif

#ifdef DEBUG_DUPLICATE
namespace {
    // Debug state
    std::map<std::string, int> hashCollisions;
    
    uint64 getBaseValue(char c) {
        switch(c) {
            case 'A': return 7;
            case 'T': return 222;
            case 'C': return 74;
            case 'G': return 31;
            default: return 13;
        }
    }

    // Helper functions now take parameters instead of accessing private members
    bool debugBloomFilter(
        Duplicate* dup,  // Still need this for seq2intvector call
        uint64* positions, 
        const std::string& sequence,
        uint32 bufNum,
        uint64 bufLenInBits,
        uint64 bufLenInBytes,
        atomic_uchar* dupBuf
    ) {
        std::vector<bool> hashResults;
        bool finalIsDup = true;

        printf("\n=== BLOOM FILTER DEBUG for sequence: %s ===\n",
               sequence.substr(0, 20).c_str());

        for(int i = 0; i < static_cast<int>(bufNum); i++) {
            uint64 pos = positions[i] % bufLenInBits;
            uint64 bytePos = pos >> 3;
            uint32 bitOffset = pos & 0x07;
            uint8 byte = (0x01) << bitOffset;

            uint8 ret = atomic_fetch_or(dupBuf + i * bufLenInBytes + bytePos, byte);
            bool thisHashSaysDup = (ret & byte) != 0;
            hashResults.push_back(thisHashSaysDup);

            printf("  Hash %d: pos=%lu, bit_was_set=%s\n",
                   i, pos, thisHashSaysDup ? "YES" : "NO");

            finalIsDup = thisHashSaysDup;  // The bug
        }

        bool correctResult = std::all_of(hashResults.begin(), hashResults.end(),
                                       [](bool b) { return b; });

        printf("  BROKEN RESULT: %s (only last hash: %s)\n",
               finalIsDup ? "DUPLICATE" : "UNIQUE",
               hashResults.back() ? "true" : "false");
        printf("  CORRECT RESULT: %s (all hashes: %s)\n",
               correctResult ? "DUPLICATE" : "UNIQUE",
               correctResult ? "all true" : "some false");

        if(finalIsDup != correctResult) {
            printf("  *** BUG CONFIRMED: Bloom filter logic is broken! ***\n");
        }

        return finalIsDup;
    }

    void debugHashFunction(
        Duplicate* dup,  // Still need this for seq2intvector call  
        const std::string& sequence,
        uint32 bufNum,
        uint64 offsetMask,
        uint64* primeArrays
    ) {
        uint64* positions = new uint64[bufNum];
        std::fill(positions, positions + bufNum, 0ULL);

        printf("\n=== HASH ANALYSIS: %s ===\n", sequence.c_str());
        printf("Position | Char | Base | Contributions per hash buffer\n");
        printf("---------|------|------|--------------------------------\n");

        for(int p = 0; p < static_cast<int>(sequence.length()); p++) {
            uint64 base = getBaseValue(sequence[p]);
            printf("%8d | '%c' | %4lu | ", p, sequence[p], base);

            for(int i = 0; i < static_cast<int>(bufNum); i++) {
                int offset = p * static_cast<int>(bufNum) + i;
                offset &= static_cast<int>(offsetMask);
                uint64 prime = primeArrays[offset];
                uint64 contribution = prime * (base + p);
                positions[i] += contribution;

                printf("H%d:+%lu ", i, contribution);
            }
            printf("\n");
        }

        std::string hashSig;
        for(int i = 0; i < static_cast<int>(bufNum); i++) {
            hashSig += std::to_string(positions[i]);
            if(i < static_cast<int>(bufNum) - 1) hashSig += ",";
        }

        hashCollisions[hashSig]++;
        printf("Final signature: [%s]\n", hashSig.c_str());

        if(hashCollisions[hashSig] > 1) {
            printf("*** HASH COLLISION #%d with this signature! ***\n",
                   hashCollisions[hashSig]);
        }

        delete[] positions;
    }

    void debugPrimeGeneration(
        uint32 bufNum,
        uint64* primeArrays
    ) {
        printf("\n=== PRIME GENERATION ANALYSIS ===\n");

        const int sampleSize = std::min(20, static_cast<int>(bufNum * PRIME_ARRAY_LEN));
        printf("Showing first %d primes:\n", sampleSize);
        printf("Index | Prime  | Gap   | Notes\n");
        printf("------|--------|-------|--------\n");

        uint64 totalGap = 0;
        int largeGaps = 0;

        for(int i = 0; i < sampleSize; i++) {
            uint64 gap = (i > 0) ? primeArrays[i] - primeArrays[i-1] : 0;

            printf("%5d | %6lu | %5lu", i, primeArrays[i], gap);

            if(gap > 1000) {
                printf(" | *** HUGE GAP! ***");
                largeGaps++;
            }
            printf("\n");

            if(i > 0) totalGap += gap;
        }

        if(sampleSize > 1) {
            printf("\nAverage gap: %.1f\n", static_cast<double>(totalGap) / (sampleSize - 1));
            printf("Large gaps (>1000): %d/%d (%.1f%%)\n",
                   largeGaps, sampleSize - 1,
                   100.0 * largeGaps / (sampleSize - 1));
        }
    }

    void debugOffsetAliasing(
        uint32 bufNum,
        uint64 offsetMask,
        uint64* primeArrays,
        int maxSeqLen = 100
    ) {
        printf("\n=== OFFSET ALIASING ANALYSIS ===\n");
        printf("OffsetMask: 0x%lx (%lu)\n", offsetMask, offsetMask);

        std::map<int, std::vector<std::pair<int, int>>> aliasMap;

        for(int p = 0; p < maxSeqLen; p++) {
            for(int i = 0; i < static_cast<int>(bufNum); i++) {
                int offset = p * static_cast<int>(bufNum) + i;
                int masked = offset & static_cast<int>(offsetMask);
                aliasMap[masked].push_back({p, i});
            }
        }

        int aliasCount = 0;
        for(const auto& pair : aliasMap) {
            if(pair.second.size() > 1) {
                aliasCount++;
                if(aliasCount <= 5) {
                    printf("Masked[%d] -> positions: ", pair.first);
                    for(const auto& pos : pair.second) {
                        printf("(%d,%d) ", pos.first, pos.second);
                    }
                    printf("(prime: %lu)\n", primeArrays[pair.first]);
                }
            }
        }
        printf("Total aliased offsets: %d\n", aliasCount);
    }

    void debugHashAvalanche(
        Duplicate* dup,  // Still need this for seq2intvector call
        uint32 bufNum
    ) {
        printf("\n=== HASH AVALANCHE TEST ===\n");
        printf("Testing how small input changes affect hash output:\n\n");

        struct TestCase {
            std::string seq1, seq2, description;
        };

        std::vector<TestCase> tests = {
            {"AAAAAAAAAA", "AAAAAAAAAB", "Last base A->B"},
            {"AAAAAAAAAA", "AAAAABAAAA", "Middle base A->B"},
            {"ATCGATCGAT", "ATCGATCGAC", "Last base T->C"},
            {"ATCGATCGAT", "GTCGATCGAT", "First base A->G"},
            {"AAAAAAAA", "AAAAAAAAA", "Length difference"}
        };

        for(const auto& test : tests) {
            uint64* hash1 = new uint64[bufNum];
            uint64* hash2 = new uint64[bufNum];
            std::fill(hash1, hash1 + bufNum, 0ULL);
            std::fill(hash2, hash2 + bufNum, 0ULL);

            dup->seq2intvector(test.seq1.c_str(), test.seq1.length(), hash1);
            dup->seq2intvector(test.seq2.c_str(), test.seq2.length(), hash2);

            printf("Test: %s\n", test.description.c_str());
            printf("  '%s' -> ", test.seq1.c_str());
            for(int i = 0; i < static_cast<int>(bufNum); i++) printf("%lu ", hash1[i]);
            printf("\n  '%s' -> ", test.seq2.c_str());
            for(int i = 0; i < static_cast<int>(bufNum); i++) printf("%lu ", hash2[i]);

            printf("\n  Differences: ");
            bool hasAvalanche = false;
            for(int i = 0; i < static_cast<int>(bufNum); i++) {
                uint64 diff = (hash1[i] > hash2[i]) ? hash1[i] - hash2[i] : hash2[i] - hash1[i];
                printf("%lu ", diff);
                if(diff > hash1[i] / 10) hasAvalanche = true;
            }
            printf("%s\n\n", hasAvalanche ? "" : "*** NO AVALANCHE! ***");

            delete[] hash1;
            delete[] hash2;
        }
    }

} // end anonymous namespace

// The friend function - this CAN access private members
void runAllDebugTests(Duplicate* dup) {
    printf("\n==============================================================\n");
    printf("DUPLICATE DETECTION DEBUG ANALYSIS\n");
    printf("==============================================================\n");

    // Extract private members once (friend function can access these)
    uint32 bufNum = dup->mBufNum;
    uint64 bufLenInBits = dup->mBufLenInBits;
    uint64 bufLenInBytes = dup->mBufLenInBytes;
    uint64 offsetMask = dup->mOffsetMask;
    atomic_uchar* dupBuf = dup->mDupBuf;
    uint64* primeArrays = dup->mPrimeArrays;

    // Pass these values to helper functions
    debugPrimeGeneration(bufNum, primeArrays);
    debugOffsetAliasing(bufNum, offsetMask, primeArrays);
    debugHashAvalanche(dup, bufNum);  // Still pass dup for seq2intvector call

    std::vector<std::string> testSeqs = {
        "ATCGATCGATCG",
        "ATCGATCGATCG", // Exact duplicate
        "ATCGATCGATCA", // One base different
        "AAAAAAAAAA",   // Homopolymer
        "ABABABABAB"    // Repetitive
    };

    for(const auto& seq : testSeqs) {
        debugHashFunction(dup, seq, bufNum, offsetMask, primeArrays);

        uint64* positions = new uint64[bufNum];
        std::fill(positions, positions + bufNum, 0ULL);
        dup->seq2intvector(seq.c_str(), seq.length(), positions);
        debugBloomFilter(dup, positions, seq, bufNum, bufLenInBits, bufLenInBytes, dupBuf);
        delete[] positions;
    }

    printf("\n==============================================================\n");
    printf("DEBUG ANALYSIS COMPLETE\n");
    printf("==============================================================\n\n");
}

#endif // DEBUG_DUPLICATE

Duplicate::Duplicate(Options* opt) {
    mOptions = opt;

    // 1G mem required
    mBufLenInBytes = 1L <<29;
    mBufNum = 2;

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
    mDupBuf = new atomic_uchar[mBufLenInBytes * mBufNum];
    if(!mDupBuf) {
        error_exit("Out of memory, failed to allocate " + to_string(mBufLenInBytes * mBufNum) + " bytes buffer for duplication analysis, please reduce the dup_accuracy_level and try again.");
    }
    memset(mDupBuf, 0, sizeof(atomic_uchar) * mBufLenInBytes * mBufNum);

    mPrimeArrays = new uint64[mBufNum * PRIME_ARRAY_LEN];
    memset(mPrimeArrays, 0, sizeof(uint64) * mBufNum * PRIME_ARRAY_LEN);
    initPrimeArrays();

    mTotalReads = 0;
    mDupReads = 0;

    #ifdef DEBUG_DUPLICATE
        runAllDebugTests(this);
    #endif
}

void Duplicate::initPrimeArrays() {
    uint64 number = 10000;
    uint64 count = 0;
    while(count < mBufNum * PRIME_ARRAY_LEN) {
        number++;
        bool isPrime = true;
        for(uint64 i=2; i<=sqrt(number); i++) {
            if(number%i == 0) {
                isPrime = false;
                break;
            }
        }
        if(isPrime) {
            mPrimeArrays[count] = number;
            count++;
            number += 10000;
        }
    }
}

Duplicate::~Duplicate(){
    delete[] mDupBuf;
    delete[] mPrimeArrays;
}

void Duplicate::seq2intvector(const char* data, int len, uint64* output, int posOffset) {
    for(int p=0; p<len; p++) {
        uint64 base = 0;
        switch(data[p]) {
            case 'A':
                base = 7;
                break;
            case 'T':
                base = 222;
                break;
            case 'C':
                base = 74;
                break;
            case 'G':
                base = 31;
                break;
            default:
                base = 13;
        }
        for(int i=0; i<mBufNum; i++) {
            int offset = (p+posOffset)*mBufNum + i;
            offset &= mOffsetMask;
            output[i] += mPrimeArrays[offset] * (base + (p+posOffset));
        }
    }
}

bool Duplicate::checkRead(Read* r) {
    uint64* positions = new uint64[mBufNum];

    // init
    for(int i=0; i<mBufNum; i++)
        positions[i] = 0;
    int len = r->length();
    seq2intvector(r->mSeq->c_str(), len, positions);
    bool isDup = applyBloomFilter(positions);
    delete[] positions;

    mTotalReads++;
    if(isDup)
        mDupReads++;

    return isDup;
}

bool Duplicate::checkPair(Read* r1, Read* r2) {
    uint64* positions = new uint64[mBufNum];
    
    // init
    for(int i=0; i<mBufNum; i++)
        positions[i] = 0;
    seq2intvector(r1->mSeq->c_str(), r1->length(), positions);
    seq2intvector(r2->mSeq->c_str(), r2->length(), positions, r1->length());
    bool isDup = applyBloomFilter(positions);
    delete[] positions;

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
        uint8 ret = atomic_fetch_or(mDupBuf + i * mBufLenInBytes + bytePos, byte);
        isDup = (ret & byte) != 0;
    }
    return isDup;
}

double Duplicate::getDupRate() {
    if(mTotalReads == 0)
        return 0.0;
    return (double)mDupReads/(double)mTotalReads;
}
