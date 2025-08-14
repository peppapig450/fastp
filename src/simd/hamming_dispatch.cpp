#include <cstddef>
#include <cstdint>

#include "bit.h"

#if defined(__AVX512BW__) || defined(__AVX2__) || defined(__SSE2__)
    #include <immintrin.h>
#endif
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    #include <arm_neon.h>
#endif

#if defined(_MSC_VER)
    #include <intrin.h>  // __cpuid, __cpuidex

    #include <array>
#endif
#if defined(__arm__) && defined(__linux__)
    #include <asm/hwcap.h>  // HWCAP_* flags
    #include <sys/auxv.h>   // getuaxval
#endif

namespace simd {
namespace detail {

// Helper macro for per-function ISA attributes on GCC/clang
#if defined(__GNUC__) || defined(__clang__)
    #define SIMD_TARGET(x) __attribute__((target(x)))
#else
    #define SIMD_TARGET(x)
#endif

namespace {

auto scalar_impl(const char* lhsBytes, const char* rhsBytes, std::size_t length, int maxMismatches)
    -> int {
    if (lhsBytes == nullptr || rhsBytes == nullptr) {
        return -1;
    }

    if (length == 0) {
        return 0;  // Early exit: there is no bytes to compare
    }

    int mismatchCount = 0;
    for (std::size_t i = 0; i < length; ++i) {
        mismatchCount += static_cast<int>(lhsBytes[i] != rhsBytes[i]);
        if (mismatchCount > maxMismatches) {
            return mismatchCount;
        }
    }

    return mismatchCount;
}

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
SIMD_TARGET("sse2")
auto sse2_impl(const char* lhsBytes, const char* rhsBytes, std::size_t length, int maxMismatches)
    -> int {
    if (lhsBytes == nullptr || rhsBytes == nullptr) {
        return -1;
    }

    if (length == 0) {
        return 0;
    }

    constexpr std::size_t   SSE2_REGISTER_BYTES = 16U;  // 128 bits / 8
    constexpr std::uint32_t MASK16              = 0xFFFFU;

    using simd::bit::popcount;
    int         mismatchCount = 0;
    std::size_t ptrOffset     = 0;

    for (; ptrOffset + SSE2_REGISTER_BYTES <= length; ptrOffset += SSE2_REGISTER_BYTES) {
        __m128i lhsVec    = _mm_loadu_si128(reinterpret_cast<const __m128i*>(lhsBytes + ptrOffset));
        __m128i rhsVec    = _mm_loadu_si128(reinterpret_cast<const __m128i*>(rhsBytes + ptrOffset));
        __m128i cmpEqMask = _mm_cmpeq_epi8(lhsVec, rhsVec);

        auto movemask = static_cast<std::uint32_t>(_mm_movemask_epi8(cmpEqMask));
        auto mask     = (movemask ^ MASK16);

        mismatchCount += popcount(mask);
        if (mismatchCount > maxMismatches) {
            return mismatchCount;
        }
    }

    // Scalar tail fallback
    for (; ptrOffset < length; ++ptrOffset) {
        mismatchCount += static_cast<int>(lhsBytes[ptrOffset] != rhsBytes[ptrOffset]);
        if (mismatchCount > maxMismatches) {
            return mismatchCount;
        }
    }

    return mismatchCount;
}
#endif

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
SIMD_TARGET("avx2")
auto avx2_impl(const char* lhsBytes, const char* rhsBytes, std::size_t length, int maxMismatches)
    -> int {
    if (lhsBytes == nullptr || rhsBytes == nullptr) {
        return -1;
    }

    if (length == 0) {
        return 0;
    }

    constexpr std::size_t AVX2_REGISTER_BYTES = 32U;  // 256 bits / 8

    using simd::bit::popcount;
    int         mismatchCount = 0;
    std::size_t ptrOffset     = 0;

    for (; ptrOffset + AVX2_REGISTER_BYTES <= length; ptrOffset += AVX2_REGISTER_BYTES) {
        __m256i lhsVec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(lhsBytes + ptrOffset));
        __m256i rhsVec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(rhsBytes + ptrOffset));
        __m256i cmpEq  = _mm256_cmpeq_epi8(lhsVec, rhsVec);

        auto cmpEqMask = static_cast<std::uint32_t>(_mm256_movemask_epi8(cmpEq));

        mismatchCount += static_cast<int>(AVX2_REGISTER_BYTES) - popcount(cmpEqMask);
        if (mismatchCount > maxMismatches) {
            return mismatchCount;
        }
    }

    // tail: reuse SSE2 or scalar
    #if defined(__GNUC__) || defined(__clang__)
    return mismatchCount
           + sse2_impl(lhsBytes + ptrOffset,
                       rhsBytes + ptrOffset,
                       length - ptrOffset,
                       maxMismatches - mismatchCount);
    #else
    return mismatchCount
           + scalar_impl(lhsBytes + ptrOffset,
                         rhsBytes + ptrOffset,
                         length - ptrOffset,
                         maxMismatches - mismatchCount);
    #endif
}
#endif

#if defined(__x86_64__) || defined(_M_X86)
SIMD_TARGET("avx512bw")
auto avx512bw_impl(const char* lhsBytes,
                   const char* rhsBytes,
                   std::size_t length,
                   int         maxMismatches) -> int {
    if (lhsBytes == nullptr || rhsBytes == nullptr) {
        return -1;
    }

    if (length == 0) {
        return 0;
    }

    constexpr std::size_t AVX512_REGISTER_BYTES = 64U;  // 512 bits / 8

    using simd::bit::popcount;
    int         mismatchCount = 0;
    std::size_t ptrOffset     = 0;

    for (; ptrOffset + AVX512_REGISTER_BYTES <= length; ptrOffset += AVX512_REGISTER_BYTES) {
        __m512i   lhsVec = _mm512_loadu_si512(lhsBytes + ptrOffset);
        __m512i   rhsVec = _mm512_loadu_si512(rhsBytes + ptrOffset);
        __mmask64 cmpEq  = _mm512_cmpeq_epi8_mask(lhsVec, rhsVec);

        // Count mismatches = 64 - #equal lanes
        mismatchCount +=
            static_cast<int>(AVX512_REGISTER_BYTES) - popcount(static_cast<std::uint64_t>(cmpEq));
        if (mismatchCount > maxMismatches) {
            return mismatchCount;
        }
    }

    // Masked tail (0..63 bytes remain)
    const std::size_t rem = length - ptrOffset;
    if (rem != 0U) {
        const __mmask64 tailMask = (1ULL << rem) - 1ULL;  // rem in [0,63], so shift is safe
        const __m512i   lhsVec   = _mm512_maskz_loadu_epi8(tailMask, lhsBytes + ptrOffset);
        const __m512i   rhsVec   = _mm512_maskz_loadu_epi8(tailMask, rhsBytes + ptrOffset);

        // Ignore zeroed lanes: AND with tailMask
        const std::uint64_t eqMask =
            static_cast<std::uint64_t>(_mm512_cmpeq_epi8_mask(lhsVec, rhsVec)) & tailMask;

        mismatchCount += static_cast<int>(rem) - popcount(eqMask);
    }

    return mismatchCount;
}
#endif

using impl_t = decltype(&scalar_impl);

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
auto choose_best_impl(bool has_avx512bw, bool has_avx2, bool has_sse2) -> impl_t {
    if (has_avx512bw) {
        return &avx512bw_impl;
    }

    if (has_avx2) {
        return &avx2_impl;
    }

    if (has_sse2) {
        return &sse2_impl;
    }

    return &scalar_impl;
}
#endif

#if defined (__arm__) || defined(__aarch64__)
    #if defined(__GNUC__) || defined (__clang__)
        #if defined (__aarch64__)
            #define NEON_ATTR // AArch64 has NEON baseline
        #else
            #define NEON_ATTR __attribute__((target("fpu=neon")))
        #endif
    #else
        #define NEON_ATTR
    #endif

NEON_ATTR
static auto neon_impl(const char* lhsBytes, const char* rhsBytes, std::size_t length, int maxMismatches) -> int {
    if (lhsBytes == nullptr || rhsBytes == nullptr) {
        return -1;
    }

    if (length == 0) {
        return 0;
    }

    constexpr std::size_t NEON_REGISTER_BYTES = 16U; // 128-bits / 8

    int mismatchCount = 0;
    std::size_t ptrOffset = 0;

    for (; ptrOffset + NEON_REGISTER_BYTES <= length; ptrOffset += NEON_REGISTER_BYTES) {
        uint8x16_t lhsVec = vld1q_u8(reinterpret_cast<const std::uint8_t*>(lhsBytes + ptrOffset));
        uint8x16_t rhsVec = vld1q_u8(reinterpret_cast<const std::uint8_t*>(rhsBytes + ptrOffset));

        // eq = 0xFF where equal, 0x00 where not
        uint8x16_t eq = vceqq_u8(lhsVec, rhsVec);

        // ne = 0xFF where mismatch, 0x00 where not
        uint8x16_t cmpNe = vmvnq_u8(eq);

        // Convert 0xFF/0x00 -> 1/0 per lane, then sum lanes
        uint8x16_t ones = vshrq_n_u8(cmpNe, 7);

        #if defined (__aarch64__)
            mismatchCount += static_cast<int>(vaddvq_u8(ones));
        #else
            // Portable reduction for 32-bit ARM
            uint16x8_t sum16 = vpaddlq_u8(ones); // pairwise add -> u16
            uint32x4_t sum32 = vpaddlq_u16(sum16); // -> u32
            uint64x2_t sum64 = vpaddlq_u32(sum32); // -> u64

            mismatchCount += static_cast<int>(vgetq_lane_u64(sum64, 0) + vgetq_lane_u64(sum64, 1));
        #endif

        if (mismatchCount > maxMismatches) {
            return mismatchCount;
        }
    }

    // Scalar tail loop
    for (; ptrOffset < length; ++ptrOffset) {
        mismatchCount += static_cast<int>(lhsBytes[ptrOffset] != rhsBytes[ptrOffset]);
        if (mismatchCount > maxMismatches) {
            return mismatchCount;
        }
    }

    return mismatchCount;
}
#endif

}  // namespace

// Bit masks for CPUID features (MSVC path)
// These are documented at the following:
// see: https://shell-storm.org/x86doc/CPUID.html
// and: https://www.felixcloutier.com/x86/cpuid
#if defined(_MSC_VER)
constexpr int EDX_SSE2_BIT     = 1U << 26U;  // leaf 1, EDX
constexpr int EBX_AVX2_BIT     = 1U << 5U;   // leaf 7, subleaf 0, EBX
constexpr int EBX_AVX512F_BIT  = 1U << 16U;  // leaf 7, subleaf 0, EBX
constexpr int EBX_AVX512BW_BIT = 1U << 30U;  // leaf 7, subleaf 0, EBX
#endif

using HammingCapFunc = int (*)(const char*, const char*, std::size_t, int);

auto pick_impl() -> HammingCapFunc {
#if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
    // GCC/clang on x86: we can use builtins with readable names
    const bool has_avx512bw = __builtin_cpu_supports("avx512bw");
    const bool has_avx2     = __builtin_cpu_supports("avx2");
    const bool has_sse2     = __builtin_cpu_supports("sse2");

    return choose_best_impl(has_avx512bw, has_avx2, has_sse2);

#elif defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
    // MSVC on x86: call CPUID and use named bit masks
    std::array<int, 4> regs_leaf1 {{0}};  // EAX, EBX, ECX, EDX
    __cpuid(regs_leaf1.data(), 1);
    const bool has_sse2 = (regs_leaf1[3] & EDX_SSE2_BIT) != 0;

    std::array<int, 4> regs_leaf7_s0 {{0}};  // leaf 7, subleaf 0
    __cpuidex(regs_leaf7_s0.data(), 7, 0);
    const bool has_avx2     = (regs_leaf7_s0[1] & EBX_AVX2_BIT) != 0;
    const bool has_avx512f  = (regs_leaf7_s0[1] & EBX_AVX512F_BIT) != 0;
    const bool has_avx512bw = (regs_leaf7_s0[1] & EBX_AVX512BW_BIT) != 0 && has_avx512f;

    return choose_best_impl(has_avx512bw, has_avx2, has_sse2);

#elif defined(__aarch64__)
    // AArch64: NEON is baseline
    return &neon_impl;

#elif defined(__arm__)
    #if defined(__linux__)
    // ARMv7 Linux: check HWCAP for NEON
    unsigned long hwcaps = getauxval(AT_HWCAP);
    if (hwcaps & HWCAP_NEON) {
        return &neon_impl;
    }

    return &scalar_impl;
    #else
    // Best effort on all other ARMv7 platforms
    return &neon_impl;
    #endif

#else
    // Fallback
    return &scalar_impl;
#endif
}
}  // namespace detail

auto hammingCap(const char* lhsBytes,
                const char* rhsBytes,
                std::size_t length,
                int         maxMismatches) noexcept -> int {
    static detail::HammingCapFunc impl = detail::pick_impl();
    return impl(lhsBytes, rhsBytes, length, maxMismatches);
}

}  // namespace simd
