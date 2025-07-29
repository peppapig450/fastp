#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

class Matcher {
public:
    Matcher()  = default;
    ~Matcher() = default;

    // Delete the copy constructor and assignment operator as they're unneeded
    Matcher(const Matcher&)                    = delete;
    auto operator=(const Matcher&) -> Matcher& = delete;

    static auto matchWithOneInsertion(const char* insertionData,
                                      const char* normalData,
                                      int         compareLength,
                                      int         diffLimit) -> bool;
    static auto diffWithOneInsertion(const char* insertionData,
                                     const char* normalData,
                                     int         compareLength,
                                     int         diffLimit) -> int;

    /**
     * Landau-Vishkin edit-distance (unit-cost Levenshtein) limited to
     * @p maxEdits.
     *
     * @return
     *   * the minimal distance ∈ [0,maxAllowedEdits] when it does not
     *     exceed the bound or
     *   * -1 when the distance > @p maxAllowedEdits.
     */
    static auto editDistance(const char* firstString,
                             int         lenA,
                             const char* secondString,
                             int         lenB,
                             int         maxAllowedEdits) -> int;

    static auto test() -> bool;

private:
    // Buffer size threshold for using stack allocation.
    // 1024 ints ~= 4 KiB per buffer (total 8 KiB), safe for typical stack usage.
    static constexpr std::size_t StackLimit = 1024;

    // Template function to handle stack vs heap allocation.
    //
    // As long as `compareLength` is less than or equal to `StackLimit`, which is the common case,
    // we allocate buffers on the stack for performance. This avoids the overhead of heap allocation
    // while remaining within typical safe stack usage.
    //
    // If `compareLength` exceeds the limit, we fall back to heap allocation to prevent
    // potential stack overflows.
    //
    // `func` is a lambda function (i.e., an inline anonymous function) that takes two `int*`
    // buffers: one for left mismatches, one for right mismatches. It encapsulates the core matching
    // logic.
    //
    // This design allows us to centralize memory handling while reusing the matching logic for
    // different public APIs (like `matchWithOneInsertion` and `diffWithOneInsertion`).
    template <typename ReturnType, typename Function>
    static auto allocateAndExecute(int length, Function&& func) -> ReturnType {
        if (length <= 0) {
            throw std::invalid_argument("length must be non-negative");
        }

        // This is safe because we already know that length is not negative
        // so we can cast to an unsigned long (size_t).
        auto len = static_cast<size_t>(length);

        if (len <= StackLimit) {
            // allocate our buffers on the stack
            std::array<int, StackLimit> leftMismatches;
            std::array<int, StackLimit> rightMismatches;
            // Only zero-initialize the active region [0, len).
            // Avoids the overhead of default-initializing the entire StackLimit-sized array.
            std::fill_n(leftMismatches.data(), len, 0);
            std::fill_n(rightMismatches.data(), len, 0);

            return func(leftMismatches.data(), rightMismatches.data());
        }

        // Fallback to using the heap for extremely long sequences
        // this handles the edge case and prevents stack buffer overflow bugs.
        std::vector<int> leftMismatches(len, 0);
        std::vector<int> rightMismatches(len, 0);

        return func(leftMismatches.data(), rightMismatches.data());
    }

    // Internal implementations that operate on preallocated mismatch buffers.
    // Called by the public methods via allocateAndExecute to abstract memory handling.
    static auto matchWithOneInsertionImpl(const char* insertionData,
                                          const char* normalData,
                                          int         compareLength,
                                          int         diffLimit,
                                          int*        leftMismatches,
                                          int*        rightMismatches) -> bool;
    static auto diffWithOneInsertionImpl(const char* insertionData,
                                         const char* normalData,
                                         int         compareLength,
                                         int         diffLimit,
                                         int*        leftMismatches,
                                         int*        rightMismatches) -> int;
};
