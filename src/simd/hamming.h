#pragma once

#include <cstddef>

namespace simd {

// Count byte mismatches between lhs[0..length] and rhs[0..length], early stopping at
// `maxMismatches`.
// Returns -1 if either pointer is null.
auto hammingCap(const char* lhsBytes,
                const char* rhsBytes,
                std::size_t length,
                int         maxMismatches) noexcept -> int;
}  // namespace simd
