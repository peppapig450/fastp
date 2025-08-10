#include <benchmark/benchmark.h>

#include <string_view>

#include "../local_tests/evaluator_access.hpp"
#include "fastqreader.h"
#include "include/benchmark_data.hpp"
#include "knownadapters.h"
#include "read.h"

namespace {  // anon

constexpr auto ReferenceGenomePath =
    std::string_view("reference_genome/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna");

}  // namespace

static void BM_CheckKnownAdapters_Current_Realistic(benchmark::State& state) {
    const auto readCount  = state.range(0);
    const auto readLength = state.range(1);

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

    for (auto _ : state) {
        auto result = EvaluatorAccess::checkKnownAdapters(evaluator, reads);
        benchmark::DoNotOptimize(result);
    }

    state.SetItemsProcessed(state.iterations() * readCount);
    state.SetBytesProcessed(state.iterations() * readCount * readLength);
}