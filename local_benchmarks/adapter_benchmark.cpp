#include <benchmark/benchmark.h>

#include <algorithm>
#include <random>
#include <string>
#include <vector>

#include "include/bench_seed.hpp"
#include "include/benchmark_utils.hpp"
#include "knownadapters.h"

// Build a synthetic corpus of sequences, optionally with an adapter at position p
static auto makeSequences(std::size_t n,
                          std::size_t len,
                          bool        include_adapter,
                          std::size_t pos      = 0,
                          double      fraction = 1.0,
                          unsigned    seed     = 880) -> std::vector<std::string> {
    std::mt19937             rng(seed);
    std::vector<std::string> seqs;
    seqs.reserve(n);

    // choose a representative adapter
    const std::string adapter = benchmark_util::pickLongestAdapter();

    std::uniform_real_distribution<double> u(0.0, 1.0);

    for (std::size_t i = 0; i < n; ++i) {
        std::string s(len, 'A');
        for (auto& b : s) {
            b = benchmark_util::randBase(rng);
        }

        if (include_adapter && u(rng) <= fraction && pos < len) {
            const auto take = std::min(len - pos, adapter.size());
            std::copy_n(adapter.data(), take, s.begin() + static_cast<long>(pos));
        }
        seqs.push_back(std::move(s));
    }
    return seqs;
}

// Exact first-match search with adapter at prefix (position 0) -> early exit expected
static void BM_AC_FindFirstAdapter_Prefix(benchmark::State& state) {
    const auto n   = static_cast<std::size_t>(state.range(0));
    const auto len = static_cast<std::size_t>(state.range(1));

    auto seqs = makeSequences(n,
                              len,
                              /*include_adapter=*/true,
                              /*pos=*/0,
                              /*fraction=*/1.0,
                              static_cast<unsigned>(bench_seed::derive_seed("BM_AC_Prefix")));
    adapters::getAhoCorasickMatcher();

    std::size_t hits = 0;
    for (auto _ : state) {
        for (const auto& s : seqs) {
            auto match = adapters::findFirstAdapter(s, 0);
            if (match.first && match.second.position == 0) {
                ++hits;
            }
        }
    }
    benchmark::DoNotOptimize(hits);
    state.SetItemsProcessed(state.iterations() * static_cast<long long>(n));
    state.SetBytesProcessed(state.iterations() * static_cast<long long>(n * len));
}

// Legacy prefix-only search using matchKnown (linear scan of adapters)
static void BM_Legacy_MatchKnown_Prefix(benchmark::State& state) {
    const auto n   = static_cast<std::size_t>(state.range(0));
    const auto len = static_cast<std::size_t>(state.range(1));

    auto seqs = makeSequences(n,
                              len,
                              /*include_adapter=*/true,
                              /*pos=*/0,
                              /*fraction=*/1.0,
                              static_cast<unsigned>(bench_seed::derive_seed("BM_Legacy_Prefix")));

    std::size_t hits   = 0;
    std::size_t misses = 0;
    for (auto _ : state) {
        for (const auto& s : seqs) {
            auto match = adapters::matchKnown(s);
            if (!match.empty()) {
                ++hits;
            } else {
                ++misses;
            }
        }
    }
    benchmark::DoNotOptimize(hits);
    benchmark::DoNotOptimize(misses);
    state.SetItemsProcessed(state.iterations() * static_cast<long long>(n));
    state.SetBytesProcessed(state.iterations() * static_cast<long long>(n * len));
}

// Aho-Corasick prefix search collecting all exact hits for comparison
static void BM_AC_MatchKnown_Prefix(benchmark::State& state) {
    const auto n   = static_cast<std::size_t>(state.range(0));
    const auto len = static_cast<std::size_t>(state.range(1));

    auto seqs =
        makeSequences(n,
                      len,
                      /*include_adapter=*/true,
                      /*pos=*/0,
                      /*fraction=*/1.0,
                      static_cast<unsigned>(bench_seed::derive_seed("BM_AC_MatchKnown_Prefix")));
    adapters::getAhoCorasickMatcher();

    std::size_t hits   = 0;
    std::size_t misses = 0;
    for (auto _ : state) {
        for (const auto& s : seqs) {
            auto matches = adapters::matchKnownAhoCorasick(s, 0);
            if (!matches.empty()) {
                ++hits;
            } else {
                ++misses;
            }
        }
    }
    benchmark::DoNotOptimize(hits);
    benchmark::DoNotOptimize(misses);
    state.SetItemsProcessed(state.iterations() * static_cast<long long>(n));
    state.SetBytesProcessed(state.iterations() * static_cast<long long>(n * len));
}

// Exact search when no adapter present -> full scan
static void BM_AC_FindFirstAdapter_NoMatch(benchmark::State& state) {
    const auto n   = static_cast<std::size_t>(state.range(0));
    const auto len = static_cast<std::size_t>(state.range(1));

    auto seqs = makeSequences(n,
                              len,
                              /*include_adapter=*/false,
                              /*pos=*/0,
                              /*fraction=*/1.0,
                              static_cast<unsigned>(bench_seed::derive_seed("BM_AC_NoMatch")));
    adapters::getAhoCorasickMatcher();

    std::size_t misses = 0;
    for (auto _ : state) {
        for (const auto& s : seqs) {
            auto match = adapters::findFirstAdapter(s, 0);
            if (!match.first) {
                ++misses;
            }
        }
    }
    benchmark::DoNotOptimize(misses);
    state.SetItemsProcessed(state.iterations() * static_cast<long long>(n));
    state.SetBytesProcessed(state.iterations() * static_cast<long long>(n * len));
}

// Multi-match (adapter at middle) using matchKnownAhoCorasick (collects all exact hits)
static void BM_AC_MatchKnown_Exact_Mid(benchmark::State& state) {
    const auto n   = static_cast<std::size_t>(state.range(0));
    const auto len = static_cast<std::size_t>(state.range(1));
    const auto pos = static_cast<std::size_t>(state.range(2));

    auto seqs = makeSequences(n,
                              len,
                              /*include_adapter=*/true,
                              pos,
                              /*fraction=*/1.0,
                              static_cast<unsigned>(bench_seed::derive_seed("BM_AC_Exact_Mid")));
    adapters::getAhoCorasickMatcher();

    std::size_t total = 0;
    for (auto _ : state) {
        for (const auto& s : seqs) {
            auto hits  = adapters::matchKnownAhoCorasick(s, 0);
            total     += hits.size();
        }
    }
    benchmark::DoNotOptimize(total);
    state.SetItemsProcessed(state.iterations() * static_cast<long long>(n));
    state.SetBytesProcessed(state.iterations() * static_cast<long long>(n * len));
}

// Approximate fallback with small mismatch budget
static void BM_AC_MatchKnown_Approximate(benchmark::State& state) {
    const auto n           = static_cast<std::size_t>(state.range(0));
    const auto len         = static_cast<std::size_t>(state.range(1));
    const auto pos         = static_cast<std::size_t>(state.range(2));
    const auto mismatches  = static_cast<int>(state.range(3));
    const auto minMatchLen = static_cast<int>(state.range(4));

    // Build reads with an adapter near-match (exact mismatches)
    // Use FASTQ generator (reuses helper that embeds adapter with mismatches)
    const auto fq =
        benchmark_util::makeNearAdapterFastq("microapprox",
                                             n,
                                             len,
                                             pos,
                                             mismatches,
                                             static_cast<unsigned>(
                                                 bench_seed::derive_seed("nearadapter")));
    auto reads = benchmark_util::loadReads(fq, n);

    // Extract just sequences to call the function directly
    std::vector<std::string> seqs;
    seqs.reserve(reads.size());
    for (auto& rptr : reads) {
        seqs.push_back(rptr->seq());
    }

    adapters::getAhoCorasickMatcher();

    std::size_t found = 0;
    for (auto _ : state) {
        for (const auto& s : seqs) {
            auto matches = adapters::matchKnownApproximate(s, mismatches, minMatchLen);
            if (!matches.empty()) {
                ++found;
            }
        }
    }
    benchmark::DoNotOptimize(found);
    state.SetItemsProcessed(state.iterations() * static_cast<long long>(n));
    state.SetBytesProcessed(state.iterations() * static_cast<long long>(n * len));
}

// clang-format off
BENCHMARK(BM_AC_FindFirstAdapter_Prefix)
    ->Args({100'000, 150})
    ->Args({100'000, 250});

BENCHMARK(BM_Legacy_MatchKnown_Prefix)
    ->Args({100'000, 150})
    ->Args({100'000, 250});

BENCHMARK(BM_AC_MatchKnown_Prefix)
    ->Args({100'000, 150})
    ->Args({100'000, 250});

BENCHMARK(BM_AC_FindFirstAdapter_NoMatch)
    ->Args({100'000, 150})
    ->Args({100'000, 250});

BENCHMARK(BM_AC_MatchKnown_Exact_Mid)
    ->Args({100'000, 150, 60})
    ->Args({100'000, 250, 100});

BENCHMARK(BM_AC_MatchKnown_Approximate)
    ->Args({10'000, 150, 40, 1, 8})
    ->Args({10'000, 150, 60, 2, 10})
    ->Args({10'000, 250, 80, 3, 12});
                                                    // clang-format on
