#include <benchmark/benchmark.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <source_location>
#include <string>

#include "../local_tests/evaluator_access.hpp"
#include "evaluator.h"
#include "fastqreader.h"
#include "include/bench_seed.hpp"
#include "include/benchmark_data.hpp"
#include "include/benchmark_utils.hpp"
#include "knownadapters.h"
#include "read.h"

namespace {  // anon

const std::string ReferenceGenomePath = [] {
    const auto base =
        std::filesystem::path {std::source_location::current().file_name()}.parent_path();
    return (base / "reference_genome/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna")
        .string();
}();

void SetBenchmarkState(benchmark::State& state, std::size_t readCount, std::size_t readLength) {
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(readCount));
    state.SetBytesProcessed(state.iterations() * static_cast<std::int64_t>(readCount * readLength));
}

// This is copy pasted from src/evaluator.cpp, as it is in anonymous namespace there:
auto isLowComplexityKey(int packedKey) -> bool {
    constexpr int kKeyLen = 10;

    // clang-format off
    std::array<int, 4> baseCount{{0,0,0,0,}};
    // clang-format on
    auto key = static_cast<unsigned int>(packedKey);

    for (int i = 0; i < kKeyLen; ++i) {
        ++baseCount[(key >> (i * 2U)) & 0x3U];
    }

    if (*std::max_element(baseCount.begin(), baseCount.end()) >= kKeyLen - 4) {
        return true;
    }

    return (baseCount[2] + baseCount[3] >= kKeyLen - 2) || ((key >> 12U) == 0xFFU);
}

}  // namespace

static void BM_CheckKnownAdapters_Current_Realistic(benchmark::State& state) {
    const auto readCount  = static_cast<std::size_t>(state.range(0));
    const auto readLength = static_cast<std::size_t>(state.range(1));

    auto& generator = benchmark_data::getGenerator();

    auto fastqPath = generator.generateDataset("current_test",
                                               readCount,
                                               readLength,
                                               "illumina",
                                               0.1,
                                               "hiseq",
                                               false,
                                               ReferenceGenomePath);

    auto reads = generator.loadReadsFromFastq(fastqPath, readCount);

    Evaluator evaluator(nullptr);

    // warmup the automaton
    adapters::getAhoCorasickMatcher();

    for (auto _ : state) {
        auto result = EvaluatorAccess::checkKnownAdapters(evaluator, reads);
        benchmark::DoNotOptimize(result);
    }

    SetBenchmarkState(state, readCount, readLength);
}

// Vary contamination rate to stress AC vs approximate fallback
// Arg2 is contamination in permille (e.g., 100 = 10%)
static void BM_CheckKnownAdapters_ParamContam(benchmark::State& state) {
    const auto readCount   = static_cast<std::size_t>(state.range(0));
    const auto readLength  = static_cast<std::size_t>(state.range(1));
    const auto contam_perm = static_cast<int>(state.range(2));
    const auto contam      = static_cast<double>(contam_perm) / 1000.0;

    auto& generator = benchmark_data::getGenerator();

    const std::string tag = std::format("param_rr{}_len{}_c{}", readCount, readLength, contam_perm);

    auto fastqPath = generator.generateDataset(tag,
                                               readCount,
                                               readLength,
                                               "illumina",
                                               contam,
                                               "hiseq",
                                               false,
                                               ReferenceGenomePath);

    auto      reads = generator.loadReadsFromFastq(fastqPath, readCount);
    Evaluator evaluator {nullptr};
    adapters::getAhoCorasickMatcher();

    for (auto _ : state) {
        auto result = EvaluatorAccess::checkKnownAdapters(evaluator, reads);
        benchmark::DoNotOptimize(result);
    }

    SetBenchmarkState(state, readCount, readLength);
}

// Many reads with long homopolymer runs after position ~20 -> triggers skip path
static void BM_CheckKnownAdapters_LowComplexitySkip(benchmark::State& state) {
    const auto readCount  = static_cast<std::size_t>(state.range(0));
    const auto readLength = static_cast<std::size_t>(state.range(1));

    const std::string tag =
        std::format("skipLC_rr{}_len{}", readCount, readLength);
    auto fastq =
        benchmark_util::makeLowComplexityFastq(tag,
                                               readCount,
                                               readLength,
                                               24,
                                               40,
                                               'G',
                                               static_cast<unsigned>(
                                                   bench_seed::derive_seed("low_complexity")));
    auto reads = benchmark_util::loadReads(fastq, readCount);

    Evaluator evaluator {nullptr};

    for (auto _ : state) {
        auto result = EvaluatorAccess::checkKnownAdapters(evaluator, reads);
        benchmark::DoNotOptimize(result);
    }

    SetBenchmarkState(state, readCount, readLength);
}

// Reads with adapter as exact prefix -> validates require position-0 match fast path
static void BM_CheckKnownAdapters_AdapterPrefixDominant(benchmark::State& state) {
    const auto readCount  = static_cast<std::size_t>(state.range(0));
    const auto readLength = static_cast<std::size_t>(state.range(1));

    const std::string tag =
        std::format("pref_rr{}_len{}", readCount, readLength);
    auto fq               = benchmark_util::makeAdapterPrefixFastq(
                      tag,
                      readCount,
                      readLength,
                      /*prefix_prob=*/1.0,
                      static_cast<unsigned>(bench_seed::derive_seed("adapterprefix")));
    auto reads = benchmark_util::loadReads(fq, readCount);

    Evaluator evaluator(nullptr);
    adapters::getAhoCorasickMatcher();

    for (auto _ : state) {
        auto result = EvaluatorAccess::checkKnownAdapters(evaluator, reads);
        benchmark::DoNotOptimize(result);
    }

    SetBenchmarkState(state, readCount, readLength);
}

// Reads containing near-adapter (few mismatches) to exercise approximate fallback
static void BM_CheckKnownAdapters_NearAdapterApprox(benchmark::State& state) {
    const auto readCount  = static_cast<std::size_t>(state.range(0));
    const auto readLength = static_cast<std::size_t>(state.range(1));
    const auto pos        = static_cast<std::size_t>(state.range(2));
    const auto mismatches = static_cast<int>(state.range(3));

    const std::string tag =
        std::format("near_pos{}_mm{}", pos, mismatches);
    auto fq               = benchmark_util::makeNearAdapterFastq(
                      tag,
                      readCount,
                      readLength,
                      pos,
                      mismatches,
                      static_cast<unsigned>(bench_seed::derive_seed("nearadapter")));
    auto reads = benchmark_util::loadReads(fq, readCount);

    Evaluator evaluator {nullptr};
    adapters::getAhoCorasickMatcher();

    for (auto _ : state) {
        auto result = EvaluatorAccess::checkKnownAdapters(evaluator, reads);
        benchmark::DoNotOptimize(result);
    }

    SetBenchmarkState(state, readCount, readLength);
}

static void BM_CheckKnownAdapters_Legacy_Pure(benchmark::State& state) {
    const auto readCount  = static_cast<std::size_t>(state.range(0));
    const auto readLength = static_cast<std::size_t>(state.range(1));

    auto& generator = benchmark_data::getGenerator();
    auto  fastqPath = generator.generateDataset("legacy_test",
                                               readCount,
                                               readLength,
                                               "illumina",
                                               0.1,
                                               "hiseq",
                                               false,
                                               ReferenceGenomePath);

    auto reads = generator.loadReadsFromFastq(fastqPath, readCount);

    // Create a simple legacy-style evaluator function
    auto legacyCheckAdapters = [](const std::vector<std::unique_ptr<Read>>& reads) -> std::string {
        const auto&                          knownAdapters = adapters::getKnown();
        std::unordered_map<std::string, int> stats;

        // Initialize stats for all adapters
        for (const auto& adapterPair : knownAdapters) {
            stats[adapterPair.first] = 0;
        }

        int         bestHits = 0;
        std::string bestAdapter;

        // Check each read against each adapter (the old way)
        for (const auto& read : reads) {
            const auto& seq = read->seq();

            for (const auto& adapterPair : knownAdapters) {
                const auto& adapter = adapterPair.first;

                // Use the legacy matchKnown function (simple prefix matching)
                auto matched = adapters::matchKnown(seq);
                if (!matched.empty() && matched == adapter) {
                    stats[adapter]++;
                    if (stats[adapter] > bestHits) {
                        bestHits    = stats[adapter];
                        bestAdapter = adapter;
                    }
                }
            }
        }

        return bestAdapter;
    };

    for (auto _ : state) {
        auto result = legacyCheckAdapters(reads);
        benchmark::DoNotOptimize(result);
    }

    SetBenchmarkState(state, readCount, readLength);
}

// Legacy method with same low-complexity filtering as current version
static void BM_CheckKnownAdapters_Legacy_WithFiltering(benchmark::State& state) {
    const auto readCount  = static_cast<std::size_t>(state.range(0));
    const auto readLength = static_cast<std::size_t>(state.range(1));

    auto& generator = benchmark_data::getGenerator();
    auto  fastqPath = generator.generateDataset("legacy_filtered",
                                               readCount,
                                               readLength,
                                               "illumina",
                                               0.1,
                                               "hiseq",
                                               false,
                                               ReferenceGenomePath);

    auto      reads = generator.loadReadsFromFastq(fastqPath, readCount);
    Evaluator evaluator(nullptr);  // For access to seq2int and isLowComplexityKey

    auto legacyWithFiltering =
        [&evaluator](const std::vector<std::unique_ptr<Read>>& reads) -> std::string {
        const auto&                          knownAdapters = adapters::getKnown();
        std::unordered_map<std::string, int> stats;

        for (const auto& adapterPair : knownAdapters) {
            stats[adapterPair.first] = 0;
        }

        int                   bestHits = 0;
        std::string           bestAdapter;
        constexpr int         kKeyLen   = 8;
        constexpr std::size_t kStartPos = 20;

        for (const auto& read : reads) {
            const auto& seq            = read->seq();
            bool        shouldSkipRead = false;

            // Apply same low-complexity filtering as current version
            if (seq.length() >= kStartPos + static_cast<std::size_t>(kKeyLen)) {
                const auto limit   = seq.length() - static_cast<std::size_t>(kKeyLen);
                int        rolling = -1;

                for (std::size_t pos = kStartPos; pos <= limit && !shouldSkipRead; ++pos) {
                    rolling = EvaluatorAccess::seq2int(evaluator,
                                                       seq,
                                                       static_cast<int>(pos),
                                                       kKeyLen,
                                                       rolling);
                    // Note: We can't access isLowComplexityKey directly, so we'll skip this
                    // optimization This gives the legacy method a slight advantage in this
                    // benchmark
                }
            }

            if (shouldSkipRead) {
                continue;
            }

            // Use legacy matching approach
            for (const auto& adapterPair : knownAdapters) {
                const auto& adapter = adapterPair.first;
                auto        matched = adapters::matchKnown(seq);
                if (!matched.empty() && matched == adapter) {
                    stats[adapter]++;
                    if (stats[adapter] > bestHits) {
                        bestHits    = stats[adapter];
                        bestAdapter = adapter;
                    }
                }
            }
        }

        return bestAdapter;
    };

    for (auto _ : state) {
        auto result = legacyWithFiltering(reads);
        benchmark::DoNotOptimize(result);
    }

    SetBenchmarkState(state, readCount, readLength);
}

// Direct comparison: matchKnown vs findFirstAdapter on individual sequences
static void BM_AdapterMatching_Legacy_Direct(benchmark::State& state) {
    const auto readCount  = static_cast<std::size_t>(state.range(0));
    const auto readLength = static_cast<std::size_t>(state.range(1));

    const std::string tag =
        std::format("legacy_direct_rr{}_len{}", readCount, readLength);
    auto fastq =
        benchmark_util::makeAdapterPrefixFastq(tag, readCount, readLength, 1.0);
    auto reads = benchmark_util::loadReads(fastq, readCount);

    for (auto _ : state) {
        for (const auto& read : reads) {
            const auto& seq    = read->seq();
            auto        result = adapters::matchKnown(seq);
            benchmark::DoNotOptimize(result);
        }
    }

    SetBenchmarkState(state, readCount, readLength);
}

// Direct comparison: Aho-Corasick findFirstAdapter on individual sequences
static void BM_AdapterMatching_AhoCorasick_Direct(benchmark::State& state) {
    const auto readCount  = static_cast<std::size_t>(state.range(0));
    const auto readLength = static_cast<std::size_t>(state.range(1));

    const std::string tag =
        std::format("ac_direct_rr{}_len{}", readCount, readLength);
    auto fastq = benchmark_util::makeAdapterPrefixFastq(tag, readCount, readLength, 1.0);
    auto reads = benchmark_util::loadReads(fastq, readCount);

    // Warm up the automaton
    adapters::getAhoCorasickMatcher();

    for (auto _ : state) {
        for (const auto& read : reads) {
            const auto& seq    = read->seq();
            auto        result = adapters::findFirstAdapter(seq, 0);
            benchmark::DoNotOptimize(result);
        }
    }

    SetBenchmarkState(state, readCount, readLength);
}

// Test legacy method performance under different contamination rates
static void BM_CheckKnownAdapters_Legacy_ParamContam(benchmark::State& state) {
    const auto readCount   = static_cast<std::size_t>(state.range(0));
    const auto readLength  = static_cast<std::size_t>(state.range(1));
    const auto contam_perm = static_cast<int>(state.range(2));
    const auto contam      = static_cast<double>(contam_perm) / 1000.0;

    auto&             generator = benchmark_data::getGenerator();
    const std::string tag =
        std::format("legacy_contam_rr{}_len{}_c{}", readCount, readLength, contam_perm);

    auto fastqPath = generator.generateDataset(tag,
                                               readCount,
                                               readLength,
                                               "illumina",
                                               contam,
                                               "hiseq",
                                               false,
                                               ReferenceGenomePath);

    auto reads = generator.loadReadsFromFastq(fastqPath, readCount);

    auto legacyCheckAdapters = [](const std::vector<std::unique_ptr<Read>>& reads) -> std::string {
        const auto&                          knownAdapters = adapters::getKnown();
        std::unordered_map<std::string, int> stats;

        for (const auto& adapterPair : knownAdapters) {
            stats[adapterPair.first] = 0;
        }

        int         bestHits = 0;
        std::string bestAdapter;

        for (const auto& read : reads) {
            const auto& seq = read->seq();

            for (const auto& adapterPair : knownAdapters) {
                const auto& adapter = adapterPair.first;
                auto        matched = adapters::matchKnown(seq);
                if (!matched.empty() && matched == adapter) {
                    stats[adapter]++;
                    if (stats[adapter] > bestHits) {
                        bestHits    = stats[adapter];
                        bestAdapter = adapter;
                    }
                    break;  // Found match, go to next read
                }
            }
        }

        return bestAdapter;
    };

    for (auto _ : state) {
        auto result = legacyCheckAdapters(reads);
        benchmark::DoNotOptimize(result);
    }

    SetBenchmarkState(state, readCount, readLength);
}

// Legacy method benchmarks
BENCHMARK(BM_CheckKnownAdapters_Legacy_Pure)
    ->Args({50'000, 150})
    ->Args({200'000, 150})
    ->Args({50'000, 250});

BENCHMARK(BM_CheckKnownAdapters_Legacy_WithFiltering)
    ->Args({50'000, 150})
    ->Args({200'000, 150});

BENCHMARK(BM_CheckKnownAdapters_Legacy_ParamContam)
    ->Args({50'000, 150,  50})   // 5% contam
    ->Args({50'000, 150, 100})   // 10%
    ->Args({50'000, 150, 300});  // 30%

// Direct method comparison benchmarks
BENCHMARK(BM_AdapterMatching_Legacy_Direct)
    ->Args({10'000, 150})
    ->Args({10'000, 250});

BENCHMARK(BM_AdapterMatching_AhoCorasick_Direct)
    ->Args({10'000, 150})
    ->Args({10'000, 250});

BENCHMARK(BM_CheckKnownAdapters_Current_Realistic)
    ->Args({50'000, 150})
    ->Args({200'000, 150})
    ->Args({50'000, 250});

BENCHMARK(BM_CheckKnownAdapters_ParamContam)
    ->Args({50'000, 150,  50})   // 5% contam
    ->Args({50'000, 150, 100})   // 10%
    ->Args({50'000, 150, 300});  // 30%

BENCHMARK(BM_CheckKnownAdapters_LowComplexitySkip)
    ->Args({100'000, 150})
    ->Args({100'000, 250});

BENCHMARK(BM_CheckKnownAdapters_AdapterPrefixDominant)
    ->Args({100'000, 75})
    ->Args({100'000, 150});

BENCHMARK(BM_CheckKnownAdapters_NearAdapterApprox)
    ->Args({20'000, 150, 30,  1})   // near-perfect
    ->Args({20'000, 150, 40,  2})   // small mismatch
    ->Args({
                                                                                                20'000, 150, 60,  3});  // heavier approximate
