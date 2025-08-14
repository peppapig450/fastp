#pragma once

#include <cstdint>
#include <type_traits>

namespace simd {
namespace bit {

template <typename IntegerType>
inline auto popcountManualFallback(IntegerType value) noexcept -> int {
    int count = 0;
    while (value != 0) {
        value &= (value - 1);
        ++count;
    }

    return count;
}

inline auto popcount32(std::uint32_t value) noexcept -> int {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcount(value);
#elif defined(_MSC_VER)
    return __popcnt(value);
#else
    // Portable fallback
    return popcountManualFallback(value);
#endif
}

inline auto popcount64(std::uint64_t value) noexcept -> int {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(value);
#elif defined(_MSC_VER) && defined(_M_X64)
    return static_cast<int>(__popcnt64(value));
#elif defined(_MSC_VER)
    // 32-bit MSVC: split in halves
    return __popcnt(static_cast<std::uint32_t>(value))
           + __popcnt(static_cast<std::uint32_t>(value >> 32U));
#else
    return popcountManualFallback(value);
#endif
}

template <class IntegerType>
inline auto popcount(IntegerType value) noexcept -> typename std::enable_if<
    std::is_integral<IntegerType>::value
        && !std::is_same<typename std::remove_cv<IntegerType>::type, bool>::value
        && (sizeof(IntegerType) <= 8),
    int>::type {
    using UnsignedInt = typename std::make_unsigned<IntegerType>::type;

    return (sizeof(UnsignedInt) <= 4U)
               ? popcount32(static_cast<std::uint32_t>(static_cast<UnsignedInt>(value)))
               : popcount64(static_cast<std::uint64_t>(static_cast<UnsignedInt>(value)));
}

}  // namespace bit
}  // namespace simd