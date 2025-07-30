#include "matcher.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <vector>

// Conditionally enable macros for compiler hints related to branch prediction
#ifndef LIKELY_BRANCH
    #if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) \
        || defined(__INTEL_LLVM_COMPILER)
        // Hint to the compiler that the condition is likely true (used for branch prediction)
        #define LIKELY_BRANCH(x) __builtin_expect(!!(x), 1)

        // Hint to the compiler that the condition is likely false (used for branch prediction)
        #define UNLIKELY_BRANCH(x) __builtin_expect(!!(x), 0)
    #else
        // MSVC and others do not support branch prediction hints so we define empty fallbacks (they
        // have no effect)
        #define LIKELY_BRANCH(x) (x)
        #define UNLIKELY_BRANCH(x) (x)
    #endif
#endif

// Helper macros for passing parameters
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

// Per compiler macros for unrolling N iterations of a loop
#ifndef UNROLL_LOOP
    #if defined(__clang__) || defined(__INTEL_LLVM_COMPILER)
        // Clang or ICX (Intel LLVM-based compiler)
        #define UNROLL_LOOP(N) _Pragma(TOSTRING(clang loop unroll_count(N)))
    #elif defined(__INTEL_COMPILER)
        // ICC classic (non-LLVM)
        #define UNROLL_LOOP(N) _Pragma(TOSTRING(unroll(N)))
    #elif defined(__GNUC__)
        #define UNROLL_LOOP(N) _Pragma(TOSTRING(GCC unroll N))
    #elif defined(_MSC_VER)
        // MSVC doesn't support unroll count reliably; fallback to simple unroll hint
        #define UNROLL_LOOP(N) __pragma(loop(unroll))
    #else
        #define UNROLL_LOOP(N)
    #endif
#endif

// Macro to apply compiler-specific restrict qualifiers for pointers.
//
// `RESTRICT` hints to the compiler that pointers annotated with it do NOT alias with other
// pointers, allowing for more aggressive optimizations (e.g., vectorization, better instruction
// scheduling).
//
// WARNING: Only apply `RESTRICT` to pointers you can *guarantee* do not alias with other pointers.
// Misusing restrict leads to undefined behavior and quiet miscompilations.
//
// TLDR: Use this on local, non-overlapping buffers (e.g., scratch space) for performance wins.
// If you’re unsure, do not use `RESTRICT`. Bad restrict is worse than no restrict.
#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) \
    || defined(__INTEL_LLVM_COMPILER)
    #define RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define RESTRICT __restrict
#else
    #define RESTRICT
#endif

static_assert(std::numeric_limits<uint64_t>::digits == 64, "This code requires 64-bit uint64_t.");

namespace {  // (anon)

constexpr std::size_t MyersLimit  = 64;  // Upper size limit for Myers' algorithm
constexpr std::size_t EqTableSize = 256;
using EqTable                     = std::array<uint64_t, EqTableSize>;

// TODO: since c++11 has no constexpr functions, generating this with a python script or something
// and then pasting might be better, as it avoids runtime cost.
/*
 * Builds a 256-entry equality (Eq) table for Myers' bit-parallel algorithm.
 * Supports patterns up to 64 bases in length.
 *
 * @param  pattern    The input pattern string.
 * @param  length     The length of the pattern. (must be <= 64)
 * @param  eqTable    Output: Eq table mapping each byte value to a bitmask.
 */
void buildEqTable(const char *pattern, int length, EqTable &eqTable) {
  assert(length >= 0 &&
         length <= MyersLimit); // Bounds check for Myers' algorithm

  std::fill(eqTable.begin(), eqTable.end(), 0ULL);

  for (std::size_t position = 0; position < static_cast<std::size_t>(length);
       ++position) {
    const auto character = static_cast<unsigned char>(pattern[position]);
    eqTable[character] |= (1ULL << position);
  }
}

}  // namespace

auto Matcher::bitapSearch64(const char *text,
                            int         textLen,
                            const char *pattern,
                            int         patLen,
                            int         diffLimit,
                            int        &posOut) -> bool {
    if (patLen == 0 || patLen > MyersLimit) {
        return false;
    }

    EqTable Eq;
    buildEqTable(pattern, patLen, Eq);

    uint64_t       VP    = ~0ULL;  // take the bitwise NOT of 0
    uint64_t       VN    = 0ULL;
    int            score = patLen;
    const uint64_t MSB   = 1ULL << (patLen - 1);  // highest bit of pattern

    for (std::size_t i = 0; i < static_cast<std::size_t>(textLen); ++i) {
        uint64_t PM = Eq[static_cast<unsigned char>(text[i])];

        // Myers bit-vector step (edit distance)
        uint64_t X  = PM | VN;
        uint64_t D0 = (((PM & VP) + VP) ^ VP) | PM;
        uint64_t HP = VN | ~(D0 | VP);
        uint64_t HN = VN & D0;

        if (HP & MSB) {
            ++score;
        } else if (HN & MSB) {
            --score;
        }

        uint64_t shiftedHP = (HP << 1) | 1ULL;
        VP                 = (HN << 1) | ~(D0 | shiftedHP);
        VN                 = D0 & shiftedHP;

        if (score <= diffLimit) {
            posOut = static_cast<int>(i) - patLen + 1;  // starting index in text
            return true;
        }
    }

    return false;
}

auto Matcher::landauVishkinSearch(const char *text,
                                  int         textLen,
                                  const char *pattern,
                                  int         patLen,
                                  int         diffLimit,
                                  int        &posOut) -> bool {
    if (patLen == 0) {
        return false;
    }

    // Banded LV for each text positiion
    const int        k = diffLimit;
    std::vector<int> V((2 * k) + 1, -1);

    // slide coarse window
    for (int i = 0; i <= textLen - patLen + k; ++i) {
        std::fill(V.begin(), V.end(), -1);
        V[k + 0] = 0;  // diagonal 0

        for (int e = 0; e <= k; ++e) {
            for (int d = -e; d <= e; d += 2) {
                int idx = k + d;
                int x;

                const int insertionCandidate = V[idx + 1];
                const int deletionCandidate  = V[idx - 1];
                if (d == -e || (d != e && deletionCandidate < insertionCandidate)) {
                    x = insertionCandidate;  // insertion
                } else {
                    x = deletionCandidate + 1;  // deletion
                }

                int y = x - d;

                // extend match on this diagonal
                while (x < patLen && y < patLen && pattern[x] == text[i + y]) {
                    ++x;
                    ++y;
                }
                V[idx] = x;

                // full pattern matched
                if (x >= patLen) {
                    posOut = i;
                    return true;
                }
            }
        }
    }
    return false;
}

auto Matcher::exactMatch(const char *text,
                         int         textLen,
                         const char *pattern,
                         int         patLen,
                         int        &posOut) noexcept -> bool {
    // nothing to do if the pattern is longer than the text
    if (patLen == 0 || patLen > textLen) {
        return false;
    }

    const int limit = textLen - patLen;
    for (int i = 0; i <= limit; ++i) {
        if (std::memcmp(text + i, pattern, patLen) == 0) {
            posOut = i;
            return true;
        }
    }

    return false;
}

// TODO: This does NOT pass the tests in adapterTrimmer
auto Matcher::locateAdapter(const char *read,
                            int         readLen,
                            const char *adapter,
                            int         adapterLen,
                            int         mismatchStride,
                            int         matchReq,
                            int        &hitPos) -> bool {
    const int maxPrefixSkip = (adapterLen >= 16)   ? 4
                              : (adapterLen >= 12) ? 3
                              : (adapterLen >= 8)  ? 2
                                                   : 0;

    for (int offset = 0; offset <= maxPrefixSkip; ++offset) {
        int patLen = adapterLen - offset;
        if (patLen < matchReq) {
            break;
        }

        const char *patPtr = adapter + offset;
        int diffLimit = patLen / mismatchStride;
        int         localPos;
        bool        adapterFound;

        // if no mismatches allowed, do a straight exact match
        if (diffLimit == 0) {
            adapterFound = exactMatch(read, readLen, patPtr, patLen, localPos);
        } else if (patLen <= MyersLimit) {
            adapterFound = bitapSearch64(read, readLen, patPtr, patLen, diffLimit, localPos);
        } else {
            adapterFound = landauVishkinSearch(read, readLen, patPtr, patLen, diffLimit, localPos);
        }

        if (adapterFound) {
            hitPos = localPos - offset;  // convert back to full-adapter origin
            return true;
        }
    }
    return false;
}

// Call the internal implementation using a lambda function.
//
// The lambda captures the parameters by value (`[=]`) and will be passed two temporary
// buffers (allocated by `allocateAndExecute`) for tracking mismatches.
// This keeps memory management separate from the algorithm logic.
auto Matcher::matchWithOneInsertion(const char *insertionData,
                                    const char *normalData,
                                    int         compareLength,
                                    int         diffLimit) -> bool {
    return allocateAndExecute<bool>(compareLength, [=](int *leftMismatches, int *rightMismatches) {
        return matchWithOneInsertionImpl(insertionData,
                                         normalData,
                                         compareLength,
                                         diffLimit,
                                         leftMismatches,
                                         rightMismatches);
    });
}

// Same pattern as above: pass a lambda to abstract allocation.
//
// The lambda captures the parameters and returns the actual difference count
// from the internal implementation.
auto Matcher::diffWithOneInsertion(const char *insertionData,
                                   const char *normalData,
                                   int         compareLength,
                                   int         diffLimit) -> int {
    return allocateAndExecute<int>(compareLength, [=](int *leftMismatches, int *rightMismatches) {
        return diffWithOneInsertionImpl(insertionData,
                                        normalData,
                                        compareLength,
                                        diffLimit,
                                        leftMismatches,
                                        rightMismatches);
    });
}

// TODO: the remaining speed improvements here likely lie with writing SIMD

// Internal implementation of matchWithOneInsertion.
//
// Computes whether the two sequences differ by at most `diffLimit` mismatches.
// *assuming exactly one insertion* has occurred in `insertionData`.
auto Matcher::matchWithOneInsertionImpl(const char   *insertionData,
                                        const char   *normalData,
                                        int           compareLength,
                                        int           diffLimit,
                                        int *RESTRICT leftMismatches,
                                        int *RESTRICT rightMismatches) -> bool {
    const auto lastIndex = compareLength - 1;

    // Step 1: Initialize the mismatch buffers at the sequence boundaries.
    // leftMismatches[0] counts mismatches at the first character.
    // rightMismatches[lastIndex] counts mismatches at the last character (excluding the inserted
    // element).
    //
    // Since `compareLength` is always > 0 (enforced by `allocateAndExecute`), these index accesses
    // are safe.
    leftMismatches[0]          = (insertionData[0] == normalData[0]) ? 0 : 1;
    rightMismatches[lastIndex] = (insertionData[compareLength] == normalData[lastIndex]) ? 0 : 1;

    const auto tailMismatch = rightMismatches[lastIndex];

    // Early exit optimization: if the mismatches at the first and last positions already exceed the
    // allowed limit, there's no point in proceeding with further comparisons.
    if (UNLIKELY_BRANCH(leftMismatches[0] + tailMismatch > diffLimit)) {
        return false;
    }

    // Step 2: Forward pass — accumulate mismatches from left to right.
    //
    // leftMismatches[i] holds the cumulative mismatch count between
    // insertionData[0..i] and normalData[0..i]. The loop breaks early if
    // at any point, the cumulative mismatches plus the tail mismatch exceed the diff limit.
    //
    // `maxValidLeft` tracks the last valid prefix boundary where the mismatch budget is not
    // exceeded.
    auto maxValidLeft = lastIndex;
    UNROLL_LOOP(8)  // Unroll this loop for 8 iterations to avoid branch mispredictions
    for (int i = 1; i < compareLength; ++i) {
        leftMismatches[i] =
            leftMismatches[i - 1] + static_cast<int>(insertionData[i] != normalData[i]);
        if (UNLIKELY_BRANCH(leftMismatches[i] + tailMismatch > diffLimit)) {
            maxValidLeft = i - 1;
            break;
        }
    }

    // Early bailout: if even the minimal prefix (first character) combined with the tail mismatch
    // exceeds the diff limit, no valid insertion point is possible — exit immediately.
    if (maxValidLeft < 0) {
        return false;
    }

    // Step 3: Backward pass — accumulate mismatches from right to left.
    //
    // rightMismatches[i] holds the cumulative mismatch count between
    // insertionData[i+1..end] and normalData[i..end-1].
    //
    // The loop terminates early when the combination of right-side mismatches and the first
    // character mismatch exceeds the allowed limit. `minValidRight` marks the earliest valid suffix
    // boundary. Instead of flooding the rightMismatches buffer with high values, we simply track
    // the valid range.
    int        minValidRight    = 0;
    bool       rightExceeded    = false;
    const auto leftEdgeMismatch = leftMismatches[0];

    // Unroll this loop for 8 iterations to avoid branch mispredictions
    UNROLL_LOOP(8)
    for (int i = lastIndex - 1; i >= 0; --i) {
        const auto nextIndex = i + 1;
        rightMismatches[i]   = rightMismatches[nextIndex]
                             + static_cast<int>(insertionData[nextIndex] != normalData[i]);
        if (UNLIKELY_BRANCH(rightMismatches[i] + leftEdgeMismatch > diffLimit)) {
            minValidRight = nextIndex;
            rightExceeded = true;
            break;
        }
    }

    // Step 4: Define the range of valid insertion points.
    //
    // Valid insertion positions are constrained by both the prefix (leftMismatches)
    // and suffix (rightMismatches) mismatch counts. Positions beyond these valid
    // ranges can be skipped altogether. If no valid range remains, we can exit early.
    const auto startPos = rightExceeded ? minValidRight : 1;
    const auto endPos   = std::min(maxValidLeft + 1, compareLength);
    if (startPos > endPos) {
        return false;
    }

    // Step 5: Evaluate feasible insertion points.
    //
    // For each candidate insertion point, we check if the combined mismatch count of the prefix
    // (leftMismatches) and suffix (rightMismatches) stays within the allowed threshold. Any
    // necessary bounds checking has been handled by previous checks, allowing direct evaluation
    // of the mismatch counts.
    //
    // The loop breaks early on a successful match without unnecessary iterations.
    for (int i = startPos; i < endPos; ++i) {
        if (LIKELY_BRANCH(leftMismatches[i - 1] + rightMismatches[i] <= diffLimit)) {
            return true;
        }
    }

    return false;
}

// Internal implementation of diffWithOneInsertion.
//
// Computes the *minimum number of mismatches* between the two sequences
// assuming exactly one insertion in `insertionData`, up to a maximum of `diffLimit`.
//
// Returns:
// - The minimum number of mismatches found (0 to diffLimit)
// - -1 if no valid insertion point meets the mismatch threshold
auto Matcher::diffWithOneInsertionImpl(const char   *insertionData,
                                       const char   *normalData,
                                       int           compareLength,
                                       int           diffLimit,
                                       int *RESTRICT leftMismatches,
                                       int *RESTRICT rightMismatches) -> int {
    const auto lastIndex = compareLength - 1;

    // Reuse the same forward and backward mismatch accumulation strategy.
    leftMismatches[0]          = (insertionData[0] == normalData[0]) ? 0 : 1;
    rightMismatches[lastIndex] = (insertionData[compareLength] == normalData[lastIndex]) ? 0 : 1;
    const auto tailMismatch    = rightMismatches[lastIndex];

    // Quick check if we can succeed at all
    if (UNLIKELY_BRANCH(leftMismatches[0] + tailMismatch > diffLimit)) {
        return -1;
    }

    // Forward pass with early termination tracking
    auto maxValidLeft = lastIndex;
    for (int i = 1; i < compareLength; ++i) {
        leftMismatches[i] = leftMismatches[i - 1] + ((insertionData[i] != normalData[i]) ? 1 : 0);

        // If extending the prefix makes mismatches exceed diffLimit,
        // then any further positions cannot be valid so stop scanning
        if (UNLIKELY_BRANCH(leftMismatches[i] + tailMismatch > diffLimit)) {
            maxValidLeft = i - 1;
            break;
        }
    }

    // Early exit if no valid positions available
    if (UNLIKELY_BRANCH(maxValidLeft < 1)) {
        return -1;
    }

    auto minValidRight = 0;
    auto rightExceeded = false;

    // Process backward only until we find invalid positions
    const auto leftMatchesEdge = leftMismatches[0];
    for (int i = lastIndex - 1; i >= 0; --i) {
        const auto nextIndex = i + 1;
        rightMismatches[i] =
            rightMismatches[nextIndex] + ((insertionData[nextIndex] != normalData[i]) ? 1 : 0);

        if (UNLIKELY_BRANCH(rightMismatches[i] + leftMatchesEdge > diffLimit)) {
            minValidRight = nextIndex;
            rightExceeded = true;
            // Instead of filling array, just mark the boundary
            break;
        }
    }

    // Find minumum difference only in valid range
    int minimumDifference = std::numeric_limits<int>::max();

    // Determine actual range to check
    const auto startPos = rightExceeded ? minValidRight : 1;
    const auto endPos   = std::min(maxValidLeft + 1, compareLength);

    // Early exit if no overlap
    if (UNLIKELY_BRANCH(startPos >= endPos)) {
        return -1;
    }

    // Single pass through valid positions only
    for (int i = startPos; i < endPos; ++i) {
        const auto totalDifferences = leftMismatches[i - 1] + rightMismatches[i];

        if (LIKELY_BRANCH(totalDifferences <= diffLimit)) {
            minimumDifference = std::min(minimumDifference, totalDifferences);

            // If we found a perfect match, no need to continue further
            if (UNLIKELY_BRANCH(minimumDifference == 0)) {
                return 0;
            }
        }
    }

    return (minimumDifference == std::numeric_limits<int>::max()) ? -1 : minimumDifference;
}

bool Matcher::test() {
    const char *normal  = "ACGTAC";
    const char *withIns = "ACGTTAC";  // insert T in the middle

    if (!matchWithOneInsertion(withIns, normal, 6, 1)) {
        std::cerr << "matchWithOneInsertion failed for insertion case" << std::endl;
        return false;
    }
    int diff = diffWithOneInsertion(withIns, normal, 6, 1);
    if (diff != 0) {
        std::cerr << "diffWithOneInsertion expected 0 but got " << diff << " for insertion case"
                  << std::endl;
        return false;
    }

    const char *withMismatch = "ACGTTAG";  // one mismatch after insertion
    if (!matchWithOneInsertion(withMismatch, normal, 6, 1)) {
        std::cerr << "matchWithOneInsertion failed for mismatch case" << std::endl;
        return false;
    }
    diff = diffWithOneInsertion(withMismatch, normal, 6, 2);
    if (diff != 1) {
        std::cerr << "diffWithOneInsertion expected 1 but got " << diff << " for mismatch case"
                  << std::endl;
        return false;
    }

    // should fail when mismatch exceeds allowed limit
    if (matchWithOneInsertion(withMismatch, normal, 6, 0)) {
        std::cerr << "matchWithOneInsertion unexpectedly succeeded when diffLimit=0" << std::endl;
        return false;
    }
    diff = diffWithOneInsertion(withMismatch, normal, 6, 0);
    if (diff != -1) {
        std::cerr << "diffWithOneInsertion expected -1 but got " << diff << " when diffLimit=0"
                  << std::endl;
        return false;
    }

    return true;
}
