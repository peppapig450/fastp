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
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(readCount * readLength));
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

    auto fastq = benchmark_util::makeLowComplexityFastq("skipLC", readCount, readLength);
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

    auto fq =
        benchmark_util::makeAdapterPrefixFastq("pref", readCount, readLength, /*prefix_prob=*/1.0);
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

    auto fq = benchmark_util::makeNearAdapterFastq("near", readCount, readLength, pos, mismatches);
    auto reads = benchmark_util::loadReads(fq, readCount);

    Evaluator evaluator {nullptr};
    adapters::getAhoCorasickMatcher();

    for (auto _ : state) {
        auto result = EvaluatorAccess::checkKnownAdapters(evaluator, reads);
        benchmark::DoNotOptimize(result);
    }

    SetBenchmarkState(state, readCount, readLength);
}

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
