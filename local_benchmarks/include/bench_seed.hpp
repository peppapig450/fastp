#pragma once

#include <array>
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
 * Creates a properly seeded Mersenne Twister generator.
 *
 * The generator is seeded using a high-entropy seed_seq constructed from
 * the derived 64-bit seed to ensure proper initialization of the MT's
 * internal state.
 *
 * @param label String identifier for reproducible seeding
 * @param stream_id Optional stream identifier for multiple generators
 * @return Fully initialized MT19937 generator
 */
[[nodiscard]] inline auto make_generator(std::string_view label, std::uint64_t stream_id = 0)
    -> std::mt19937 {
    const auto seed = derive_seed(label, stream_id);

    // Split 64-bit seed into 32-bit words for proper seed_seq initialization
    const std::array<std::uint32_t, 2> seed_words {
        {static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32UL)}};

    std::seed_seq sequence(seed_words.begin(), seed_words.end());
    return std::mt19937 {sequence};
}

/**
 * Creates a generator for an explicit numeric seed (primarily for testing)
 *
 * @param explicit_seed Direct 64-bit seed value
 * @return Seeded MT19937 generator
 */
[[nodiscard]] inline auto make_generator(std::uint64_t explicit_seed) -> std::mt19937 {
    const std::array<std::uint32_t, 2> seed_words {
        {static_cast<std::uint32_t>(explicit_seed),
         static_cast<std::uint32_t>(explicit_seed >> 32UL)}};

    std::seed_seq sequence(seed_words.begin(), seed_words.end());
    return std::mt19937 {sequence};
}

}  // namespace bench_seed
