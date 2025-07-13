#include "matcher.h"

#include <algorithm>
#include <iostream>
#include <limits>

// Conditionally enable macros for compiler hints related to branch prediction
#ifndef LIKELY_BRANCH
    #if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
        // Hint to the compiler that the condition is likely true (used for branch prediction)
        #define LIKELY_BRANCH(x) __builtin_expect(!!(x), 1)

        // Hint to the compiler that the condition is likely false (used for branch prediction)
        #define UNLIKELY_BRANCH(x) __builtin_expect(!!(x), 0)
    #else
        // MSVC and others do not support branch prediction hints so we define empty fallbacks (they have no effect)
        #define LIKELY_BRANCH(x) (x)
        #define UNLIKELY_BRANCH(x) (x)
    #endif
#endif

// Call the internal implementation using a lambda function.
//
// The lambda captures the parameters by value (`[=]`) and will be passed two temporary
// buffers (allocated by `allocateAndExecute`) for tracking mismatches.
// This keeps memory management separate from the algorithm logic.
auto Matcher::matchWithOneInsertion(const char *insertionData, const char *normalData, int compareLength, int diffLimit)
    -> bool {
    return allocateAndExecute<bool>(compareLength, [=](int *leftMismatches, int *rightMismatches) {
        return matchWithOneInsertionImpl(
            insertionData, normalData, compareLength, diffLimit, leftMismatches, rightMismatches);
    });
}

// Same pattern as above: pass a lambda to abstract allocation.
//
// The lambda captures the parameters and returns the actual difference count
// from the internal implementation.
auto Matcher::diffWithOneInsertion(const char *insertionData, const char *normalData, int compareLength, int diffLimit)
    -> int {
    return allocateAndExecute<int>(compareLength, [=](int *leftMismatches, int *rightMismatches) {
        return diffWithOneInsertionImpl(
            insertionData, normalData, compareLength, diffLimit, leftMismatches, rightMismatches);
    });
}

// TODO: the loops here can be vectorized, either with IVDEP auto-vectorization or SIMD intrinsics

// Internal implementation of matchWithOneInsertion.
//
// Computes whether the two sequences differ by at most `diffLimit` mismatches.
// *assuming exactly one insertion* has occurred in `insertionData`.
auto Matcher::matchWithOneInsertionImpl(const char *insertionData,
                                        const char *normalData,
                                        int compareLength,
                                        int diffLimit,
                                        int *leftMismatches,
                                        int *rightMismatches) -> bool {
    // Step 1: Initialize the mismatch buffers at the sequence boundaries.
    // leftMismatches[0] checks the first characters; rightMismatches[end] checks the last ones.
    //
    // The access insertionData[compareLength] is safe: it points to the extra character assumed
    // to be inserted in `insertionData`, making its length `compareLength + 1`.
    leftMismatches[0] = (insertionData[0] == normalData[0]) ? 0 : 1;

    // This index access is safe because `allocateAndExecute` throws if `compareLength` is equal to or less than 0
    rightMismatches[compareLength - 1] = (insertionData[compareLength] == normalData[compareLength - 1]) ? 0 : 1;

    // Step 2: Forward pass — accumulate mismatches from left to right.
    //
    // leftMismatches[i] contains the total number of mismatches between
    // insertionData[0..i] and normalData[0..i].
    for (int i = 1; i < compareLength; ++i) {
        leftMismatches[i] = leftMismatches[i - 1];
        if (insertionData[i] != normalData[i]) {
            ++leftMismatches[i];
        }

        // Early exit: if the total mismatches (so far + final tail mismatch)
        // exceed the allowed limit, we can skip further checking.
        if (leftMismatches[i] + rightMismatches[compareLength - 1] > diffLimit) {
            break;
        }
    }

    // Step 3: Backward pass — accumulate mismatches from right to left.
    //
    // rightMismatches[i] stores mismatches from insertionData[i+1..end] and normalData[i..end-1].
    for (int i = compareLength - 2; i >= 0; --i) {
        rightMismatches[i] = rightMismatches[i + 1];
        if (insertionData[i + 1] != normalData[i]) {
            ++rightMismatches[i];
        }

        // If we exceed the diff limit during the backward pass, we propagate a high value
        // to prevent false positives during insertion point evaluation.
        if (rightMismatches[i] + leftMismatches[0] > diffLimit) {
            for (int rmIdx = 0; rmIdx < i; ++rmIdx) {
                rightMismatches[rmIdx] = diffLimit + 1;
            }
            break;
        }
    }

    // Inside the loop this is repeatedly calculated despite not changing, so we hoist it
    // outside for optimization.
    const auto tailMismatch = rightMismatches[compareLength - 1];

    // Check if we can possibly succeed before we enter the loop
    if (leftMismatches[0] + tailMismatch > diffLimit) {
        return false;
    }

    // Step 4: Evaluate all valid insertion points.
    //
    // At each index `i`, we check whether the prefix (up to i-1) and the suffix (from i)
    // match within the allowed mismatch budget. If so, return true.
    for (int i = 1; i < compareLength; ++i) {
        const auto totalDifferences = leftMismatches[i - 1] + rightMismatches[i];
        if (totalDifferences <= diffLimit) {
            return true;
        }

        // Check if subsequent iterations can possibly succeed
        if (i < compareLength - 1 && leftMismatches[i] + tailMismatch > diffLimit) {
            return false;
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
auto Matcher::diffWithOneInsertionImpl(const char *insertionData,
                                       const char *normalData,
                                       int compareLength,
                                       int diffLimit,
                                       int *leftMismatches,
                                       int *rightMismatches) -> int {
    const auto lastIndex = compareLength - 1;

    // Reuse the same forward and backward mismatch accumulation strategy.
    leftMismatches[0] = (insertionData[0] == normalData[0]) ? 0 : 1;
    rightMismatches[lastIndex] = (insertionData[compareLength] == normalData[lastIndex]) ? 0 : 1;
    const auto tailMismatch = rightMismatches[lastIndex];

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
    for (int i = lastIndex - 1; i >= 0; --i) {
        rightMismatches[i] = rightMismatches[i + 1] + ((insertionData[i + 1] != normalData[i]) ? 1 : 0);

        if (UNLIKELY_BRANCH(rightMismatches[i] + leftMismatches[0] > diffLimit)) {
            minValidRight = i + 1;
            rightExceeded = true;
            // Instead of filling array, just mark the boundary
            break;
        }
    }

    // Find minumum difference only in valid range
    int minimumDifference = std::numeric_limits<int>::max();

    // Determine actual range to check
    const auto startPos = rightExceeded ? minValidRight : 1;
    const auto endPos = std::min(maxValidLeft + 1, compareLength);

    // Early exit if no overlap
    if (UNLIKELY_BRANCH(startPos >= endPos)) {
        return -1;
    }

    // Single pass through valid positions only
    for (int i = startPos; i < endPos; ++i) {
        const auto isRightInvalid = rightExceeded && i < minValidRight;

        // Skip positions we know are invalid from backward pass
        if (UNLIKELY_BRANCH(isRightInvalid)) {
            continue;
        }

        // For positions before minValidRight when rightExceeded=true,
        // we know that rightMismatches[i] would be greater than diffLimit.
        const auto rightVal = (isRightInvalid ? (diffLimit + 1) : rightMismatches[i]);

        const auto totalDifferences = leftMismatches[i - 1] + rightVal;

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
    const char *normal = "ACGTAC";
    const char *withIns = "ACGTTAC";  // insert T in the middle

    if (!matchWithOneInsertion(withIns, normal, 6, 1)) {
        std::cerr << "matchWithOneInsertion failed for insertion case" << std::endl;
        return false;
    }
    int diff = diffWithOneInsertion(withIns, normal, 6, 1);
    if (diff != 0) {
        std::cerr << "diffWithOneInsertion expected 0 but got " << diff << " for insertion case" << std::endl;
        return false;
    }

    const char *withMismatch = "ACGTTAG";  // one mismatch after insertion
    if (!matchWithOneInsertion(withMismatch, normal, 6, 1)) {
        std::cerr << "matchWithOneInsertion failed for mismatch case" << std::endl;
        return false;
    }
    diff = diffWithOneInsertion(withMismatch, normal, 6, 2);
    if (diff != 1) {
        std::cerr << "diffWithOneInsertion expected 1 but got " << diff << " for mismatch case" << std::endl;
        return false;
    }

    // should fail when mismatch exceeds allowed limit
    if (matchWithOneInsertion(withMismatch, normal, 6, 0)) {
        std::cerr << "matchWithOneInsertion unexpectedly succeeded when diffLimit=0" << std::endl;
        return false;
    }
    diff = diffWithOneInsertion(withMismatch, normal, 6, 0);
    if (diff != -1) {
        std::cerr << "diffWithOneInsertion expected -1 but got " << diff << " when diffLimit=0" << std::endl;
        return false;
    }

    return true;
}
