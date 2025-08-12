#pragma once

#include <array>
#include <concepts>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <ranges>
#include <string_view>

namespace bench_seed {

constexpr std::uint64_t DEFAULT_SEED = 123456789ULL;

// SplitMix64 algorithm constants
constexpr std::uint64_t SPLITMIX_INCREMENT    = 0x9e3779b97f4a7c15ULL;
constexpr std::uint64_t SPLITMIX_MULTIPLIER_1 = 0xbf58476d1ce4e5b9ULL;
constexpr std::uint64_t SPLITMIX_MULTIPLIER_2 = 0x94d049bb133111ebULL;

// Hash derivation salt
constexpr std::uint64_t DERIVATION_SALT = 0x6a09e667f3bcc909ULL;

/**
 * Retrieves the base seed for benchmark random number generation.
 *
 * The seed is determined by:
 * 1. BENCH_SEED environment variable (if set)
 * 2. DEFAULT_SEED compile-time constant (fallback)
 *
 * @returns The base seed value
 */
[[nodiscard]] inline auto get_base_seed() noexcept -> std::uint64_t {
    if (const char* env_seed = std::getenv("BENCH_SEED")) {
        return std::strtoull(env_seed, nullptr, 10);
    }

    return DEFAULT_SEED;
}

/**
 * SplitMix64 pseudo-random number generator.
 *
 * A fast, high-quality PRNG used for deriving independent sub-seeds.
 * This ensures reproducible but statistically independent random streams.
 *
 * @param state Current state value
 * @return Next pseudo-random value
 */
[[nodiscard]] constexpr auto splitmix64(std::uint64_t state) noexcept -> std::uint64_t {
    state += SPLITMIX_INCREMENT;
    state  = (state ^ (state >> 30ULL)) * SPLITMIX_MULTIPLIER_1;
    state  = (state ^ (state >> 27ULL)) * SPLITMIX_MULTIPLIER_2;
    return state ^ (state >> 31ULL);
}

/**
 * Derives a deterministic 64-bit seed from a string label.
 *
 * This function creates reproducible seeds that are:
 * - Deterministic: same label + stream_id always produces same seed
 * - Independent: different labels produce statistically independent seeds
 * - Collision-resistant: small changes in input produce large changes in output
 *
 * @param label String identifier for this seed derivation
 * @param stream_id Optional stream identifier for multiple seeds from same label
 * @return Derived 64-bit seed
 */
[[nodiscard]] inline auto derive_seed(std::string_view label, std::uint64_t stream_id = 0) noexcept
    -> std::uint64_t {
    std::uint64_t hash = get_base_seed() ^ DERIVATION_SALT ^ stream_id;

    // Hash each byte of the label using SplitMix64
    for (auto byte : label | std::views::transform([](char c) noexcept {
                         return static_cast<unsigned char>(c);
                     })) {
        hash = splitmix64(hash ^ byte);
    }

    return hash;
}

/**
 * Fold/mix a 64-bit seed down to a smaller unsigned type using xor-folding.
 * This keeps good bit dispersion without call-site casts.
 */
template <std::unsigned_integral UInt>
[[nodiscard]] constexpr auto fold_seed_to(std::uint64_t source) noexcept -> UInt {
    if constexpr (sizeof(UInt) >= sizeof(std::uint64_t)) {
        // Zero-extend if the target is 64-bit or wider
        return static_cast<UInt>(source);
    } else if constexpr (sizeof(UInt) == 4) {
        auto folded  = splitmix64(source);
        folded      ^= (folded >> 32ULL);
        return static_cast<UInt>(folded);
    } else if constexpr (sizeof(UInt) == 2) {
        auto folded  = splitmix64(source);
        folded      ^= (folded >> 32ULL);
        folded      ^= (folded >> 16ULL);
        return static_cast<UInt>(folded);
    } else {  // sizeof(Uint) == 1
        auto folded  = splitmix64(source);
        folded      ^= (folded >> 32ULL);
        folded      ^= (folded >> 16ULL);
        folded      ^= (folded >> 8ULL);
        return static_cast<UInt>(folded);
    }
}

/**
 * Derive a seed directly as the unsigned type you want.
 * Example: auto s32 = derive_seed_as<std::uint32_t>("foo");
 */
template <std::unsigned_integral UInt = std::uint64_t>
[[nodiscard]] inline auto derive_seed_as(std::string_view label,
                                         std::uint64_t    stream_id = 0) noexcept -> UInt {
    return fold_seed_to<UInt>(derive_seed(label, stream_id));
}

/**
 * Helper: expand a 64-bit seed into multiple 32-bit words via SplitMix64.
 * Supplying more than two words gives seed_seq more entropy to spread
 * across engine state (works for both 32- and 64-bit Mersenne Twister variants).
 */
[[nodiscard]] inline auto expand_seed_to_u32_words(std::uint64_t base_seed) noexcept
    -> std::array<std::uint32_t, 4> {
    std::array<std::uint32_t, 4> words {};
    std::uint64_t                state = base_seed;

    // Two SplitMix64 outputs => 4x32-bit words
    state    = splitmix64(state);
    words[0] = static_cast<std::uint32_t>(state);
    words[1] = static_cast<std::uint32_t>(state >> 32ULL);

    state    = splitmix64(state);
    words[2] = static_cast<std::uint32_t>(state);
    words[3] = static_cast<std::uint32_t>(state >> 32ULL);

    return words;
}

/**
 * Creates a properly seeded random engine (templated).
 *
 * Default is std::mt19937, but you can pick std::mt19937_64 or any
 * std::uniform_random_bit_generator-compatible engine.
 *
 * Example:
 *   auto rng32 = make_generator<>("bench", 0);                // std::mt19937
 *   auto rng64 = make_generator<std::mt19937_64>("bench", 0); // 64-bit MT
 */
template <std::uniform_random_bit_generator Engine = std::mt19937>
[[nodiscard]] inline auto make_generator(std::string_view label, std::uint64_t stream_id = 0)
    -> Engine {
    const auto seed64 = derive_seed(label, stream_id);
    const auto words  = expand_seed_to_u32_words(seed64);

    std::seed_seq seq {words.begin(), words.end()};
    return Engine {seq};
}

/**
 * Creates a generator from an explicit numeric seed (templated).
 */
template <std::uniform_random_bit_generator Engine = std::mt19937>
[[nodiscard]] inline auto make_generator(std::uint64_t explicit_seed) -> Engine {
    const auto words = expand_seed_to_u32_words(explicit_seed);

    std::seed_seq seq {words.begin(), words.end()};
    return Engine {seq};
}

// These helpers are only necessary because mason2 only accepts signed 32-bit integers as seeds
// this is rather quirky, so we just mask off the sign bit so we can avoid rewriting the template
// that casts to an unsigned type.
// Upstream issue: https://github.com/seqan/seqan/issues/2570
[[nodiscard]] inline auto derive_seed_int32_nonneg(std::string_view label,
                                                   std::uint64_t    stream_id = 0) noexcept
    -> std::int32_t {
    const auto u32 = derive_seed_as<std::uint32_t>(label, stream_id);
    // Mask off sign bit so result ∈ [0, INT_MAX]; deterministic + reproducible.
    return static_cast<std::int32_t>(u32 & 0x7fffffffU);
}

// Convenience: make a Mason-friendly seed from an explicit 64-bit seed.
[[nodiscard]] inline auto fold_explicit_seed_to_int32(std::uint64_t seed64,
                                                      bool          non_negative = true) noexcept
    -> std::int32_t {
    const auto u32 = fold_seed_to<std::uint32_t>(seed64);
    return static_cast<std::int32_t>(non_negative ? (u32 & 0x7fffffffU) : u32);
}

}  // namespace bench_seed
