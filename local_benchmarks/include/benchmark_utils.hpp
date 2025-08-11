#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstring>
#include <expected>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <memory>
#include <numeric>
#include <print>
#include <random>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "benchmark_data.hpp"
#include "knownadapters.h"

class Read;  // Forward declaration

namespace fs = std::filesystem;

namespace benchmark_util {

template <class MutatorFunction>
concept Mutator = requires(MutatorFunction mutatorFunc,
                           std::string&    sequence,
                           std::string&    quality,
                           std::mt19937&   randomGenerator) {
    {mutatorFunc(sequence, quality, randomGenerator)}->std::same_as<void>;
};

struct NoopMutator {
    void operator()(std::string&, std::string&, std::mt19937&) const noexcept {}
};

inline auto workDir() -> fs::path {
    auto workingDirectory = fs::temp_directory_path() / "fastp_benchmark_util";
    fs::create_directories(workingDirectory);

    return workingDirectory;
}

inline auto mutateBase(char nucleotideBase) -> char {
    switch (nucleotideBase) {
            // clang-format off
        case 'A': return 'C';
        case 'C': return 'G';
        case 'G': return 'T';
        default:  return 'A';
            // clang-format on
    }
}

inline auto randBase(std::mt19937& rng) -> char {
    static constexpr std::array<char, 4> bases {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<int>   randomDistribution(0, 3);
    return bases[randomDistribution(rng)];
}

inline auto randQual(std::size_t sequenceLength, std::mt19937& randomGenerator) -> std::string {
    std::uniform_int_distribution<int> qualityScoreDistribution {33, 73};  // PHRED 33

    std::string qualityValues;
    qualityValues.resize_and_overwrite(sequenceLength,
                                       [&](char* bufferPtr, std::size_t bufferSize) {
                                           for (std::size_t charIdx = 0; charIdx < bufferSize;
                                                ++charIdx) {
                                               bufferPtr[charIdx] = static_cast<char>(
                                                   qualityScoreDistribution(randomGenerator));
                                           }
                                           return bufferSize;
                                       });

    return qualityValues;
}

inline auto randSeq(std::size_t sequenceLength, std::mt19937& randomGenerator) -> std::string {
    static constexpr std::array<char, 4>       nucleotideBases {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<std::size_t> baseDistribution(0, nucleotideBases.size() - 1);

    std::string sequence;
    sequence.resize_and_overwrite(sequenceLength, [&](char* bufferPtr, std::size_t bufferSize) {
        for (std::size_t charIdx = 0; charIdx < bufferSize; ++charIdx) {
            bufferPtr[charIdx] = nucleotideBases[baseDistribution(randomGenerator)];
        }
        return bufferSize;
    });

    return sequence;
}

inline auto pickLongestAdapter() -> std::string {
    const auto& adapters = adapters::getKnown();
    auto        longestAdapterIterator =
        std::ranges::max_element(adapters, {}, [](auto const& adapterPair) {
            return adapterPair.first.size();
        });

    if (longestAdapterIterator == std::ranges::end(adapters)) {
        throw std::runtime_error("No known adapters available");
    }

    return longestAdapterIterator->first;  // copy key string
}

inline auto pickMedianAdapter() -> std::string {
    const auto& adapters = adapters::getKnown();
    if (std::ranges::empty(adapters)) {
        throw std::runtime_error("No known adapters available");
    }

    std::vector<std::string_view> adapterKeys;
    adapterKeys.reserve(std::ranges::size(adapters));

    for (const auto& [key, value] : adapters) {
        adapterKeys.push_back(key);
    }

    // We take a signed midpoint and an iterator advanced by that signed difference
    const auto midpoint       = std::ranges::ssize(adapterKeys) / 2;
    const auto medianIterator = std::next(adapterKeys.begin(), midpoint);

    std::ranges::nth_element(adapterKeys, medianIterator, {}, [](std::string_view adapterSequence) {
        return adapterSequence.size();
    });

    return std::string(adapterKeys[static_cast<std::size_t>(midpoint)]);
}

template <Mutator MutatorFunction = NoopMutator>
[[nodiscard]] inline auto writeFastqWithMutatorTry(std::string_view filenameTag,
                                                   std::size_t      numReads,
                                                   std::size_t      readLength,
                                                   MutatorFunction  sequenceMutator = {},
                                                   unsigned         randomSeed      = 69)
    -> std::expected<fs::path, std::string> {
    auto outputPath = workDir() / (std::string(filenameTag) + ".fastq");
    if (fs::exists(outputPath)) {
        return outputPath;
    }

    std::ofstream output(outputPath, std::ios::binary);
    if (!output) {
        return std::unexpected("Cannot open for write: " + outputPath.string());
    }

    std::mt19937 randomGenerator(randomSeed);

    std::string sequence;
    std::string quality;
    for (std::size_t readIdx = 0; readIdx < numReads; ++readIdx) {
        sequence = randSeq(readLength, randomGenerator);
        quality  = randQual(readLength, randomGenerator);

        sequenceMutator(sequence, quality, randomGenerator);
        if (quality.size() != sequence.size()) {
            // Keep FASTQ invariant: quality length == sequence length
            quality = randQual(sequence.size(), randomGenerator);
        }

        std::print(output, "@{}_{}\n{}\n+\n{}\n", filenameTag, readIdx, sequence, quality);

        if (!output) {
            return std::unexpected("I/O error while writing: " + outputPath.string());
        }
    }

    return outputPath;
}

template <Mutator MutatorFunction = NoopMutator>
[[nodiscard]] inline auto writeFastqWithMutator(std::string_view filenameTag,
                                                std::size_t      numReads,
                                                std::size_t      readLength,
                                                MutatorFunction  sequenceMutator = {},
                                                unsigned         randomSeed      = 69) -> fs::path {
    if (auto result = writeFastqWithMutatorTry(filenameTag,
                                               numReads,
                                               readLength,
                                               std::move(sequenceMutator),
                                               randomSeed);
        result) {
        return *result;
    } else {
        throw std::runtime_error(result.error());
    }
}

// Long low-complexity (homopolymer) run starting at >= pos 20 to trigger skip path
[[nodiscard]] inline auto makeLowComplexityFastqTry(std::string_view filenameTag,
                                                    std::size_t      numReads,
                                                    std::size_t      readLength,
                                                    std::size_t      homopolymerRunStart  = 24,
                                                    std::size_t      homopolymerRunLength = 40,
                                                    char             homopolymerBase      = 'G')
    -> std::expected<fs::path, std::string> {
    return writeFastqWithMutatorTry(std::string("lowcomplex_") + std::string(filenameTag),
                                    numReads,
                                    readLength,
                                    [=](std::string& sequence,
                                        std::string& /*quality*/,
                                        std::mt19937& /*randomGenerator*/) {
                                        if (readLength > homopolymerRunStart) {
                                            const auto actualRunLength =
                                                std::min(homopolymerRunLength,
                                                         readLength - homopolymerRunStart);

                                            auto* sequenceStartPtr =
                                                sequence.data() + homopolymerRunStart;

                                            std::fill_n(sequenceStartPtr,
                                                        actualRunLength,
                                                        homopolymerBase);
                                        }
                                    });
}

[[nodiscard]] inline auto makeLowComplexityFastq(std::string_view filenameTag,
                                                 std::size_t      numReads,
                                                 std::size_t      readLength,
                                                 std::size_t      homopolymerRunStart  = 24,
                                                 std::size_t      homopolymerRunLength = 40,
                                                 char homopolymerBase = 'G') -> fs::path {
    if (auto result = makeLowComplexityFastqTry(filenameTag,
                                                numReads,
                                                readLength,
                                                homopolymerRunStart,
                                                homopolymerRunLength,
                                                homopolymerBase);
        result) {
        return *result;
    } else {
        throw std::runtime_error("make_low_complexity_fastq: " + result.error());
    }
}

[[nodiscard]] inline auto makeAdapterPrefixFastqTry(std::string_view filenameTag,
                                                    std::size_t      numReads,
                                                    std::size_t      readLength,
                                                    double           adapterPrefixProbability = 1.0,
                                                    unsigned         randomSeed               = 420)
    -> std::expected<fs::path, std::string> {
    auto adapterSequence = pickLongestAdapter();

    // Clamp to [0, 1]
    adapterPrefixProbability = std::clamp(adapterPrefixProbability, 0.0, 1.0);

    auto adapterPrefixMutator =
        [bernoulliDist = std::bernoulli_distribution(adapterPrefixProbability),
         adapterSequence =
             std::move(adapterSequence)](std::string&                  sequence,
                                         [[maybe_unused]] std::string& quality,
                                         std::mt19937& randomGenerator) mutable -> void {
        if (bernoulliDist(randomGenerator)) {
            const auto charactersToCopy = std::min(sequence.size(), adapterSequence.size());
            // Both strings are contiguous so we can use memcpy here
            std::memcpy(sequence.data(), adapterSequence.data(), charactersToCopy);
        }
    };

    return writeFastqWithMutatorTry(std::string("adapter_prefix_") + std::string(filenameTag),
                                    numReads,
                                    readLength,
                                    std::move(adapterPrefixMutator),
                                    randomSeed);
}

[[nodiscard]] inline auto makeAdapterPrefixFastq(std::string_view filenameTag,
                                                 std::size_t      numReads,
                                                 std::size_t      readLength,
                                                 double           adapterPrefixProbability = 1.0,
                                                 unsigned         randomSeed = 420) -> fs::path {
    if (auto result = makeAdapterPrefixFastqTry(filenameTag,
                                                numReads,
                                                readLength,
                                                adapterPrefixProbability,
                                                randomSeed);
        result) {
        return *result;
    } else {
        throw std::runtime_error("make_adapter_prefix_fastq: " + result.error());
    }
}

// Near-adapter with controlled mismatches at an arbitrary position
[[nodiscard]] inline auto makeNearAdapterFastqTry(std::string_view filenameTag,
                                                  std::size_t      numReads,
                                                  std::size_t      readLength,
                                                  std::size_t      adapterPosition,
                                                  int              numMismatches,
                                                  unsigned         randomSeed = 69420)
    -> std::expected<fs::path, std::string> {
    auto adapterSequence = pickMedianAdapter();

    auto nearAdapterMutator = [adapterSequence = std::move(adapterSequence),
                               adapterPosition,
                               numMismatches](std::string&                  sequence,
                                              [[maybe_unused]] std::string& quality,
                                              std::mt19937& randomGenerator) -> void {
        // If the adapter position is outside the sequence, do nothing
        if (adapterPosition >= sequence.size()) {
            return;
        }

        // Determine how many characters we can actually copy without overflowing
        const auto availableSpace   = sequence.size() - adapterPosition;
        const auto charactersToCopy = std::min<std::size_t>(availableSpace, adapterSequence.size());

        // Insert the adapter at the given position
        std::memcpy(sequence.data() + adapterPosition, adapterSequence.data(), charactersToCopy);

        std::vector<std::size_t> availableIndices;
        availableIndices.reserve(charactersToCopy);
        std::ranges::iota(availableIndices, 0);

        // Introduce exactly `numMismatches` flips (bounded by `charactersToCopy`)
        const std::size_t actualMismatches =
            std::min<std::size_t>(charactersToCopy,
                                  static_cast<std::size_t>(numMismatches > 0 ? numMismatches : 0));
        if (actualMismatches > 0) {
            std::vector<std::size_t> mutationIndices;
            mutationIndices.reserve(actualMismatches);

            std::ranges::sample(availableIndices,
                                std::back_inserter(mutationIndices),
                                static_cast<std::ptrdiff_t>(actualMismatches),
                                randomGenerator);

            for (auto mutationIndex : mutationIndices) {
                auto& nucleotideBase = sequence[adapterPosition + mutationIndex];
                nucleotideBase       = mutateBase(nucleotideBase);
            }
        }
    };

    return writeFastqWithMutatorTry(std::string("nearadapter_") + std::string(filenameTag),
                                    numReads,
                                    readLength,
                                    std::move(nearAdapterMutator),
                                    randomSeed);
}

[[nodiscard]] inline auto makeNearAdapterFastq(std::string_view filenameTag,
                                               std::size_t      numReads,
                                               std::size_t      readLength,
                                               std::size_t      adapterPosition,
                                               int              numMismatches,
                                               unsigned         randomSeed = 69420) -> fs::path {
    if (auto result = makeNearAdapterFastqTry(filenameTag,
                                              numReads,
                                              readLength,
                                              adapterPosition,
                                              numMismatches,
                                              randomSeed);
        result) {
        return *result;
    } else {
        throw std::runtime_error("make_near_adapter_fastq: " + result.error());
    }
}

inline auto loadReads(const fs::path& fastqFilePath, std::size_t maxReadsToLoad = 0)
    -> std::vector<std::unique_ptr<Read>> {
    return benchmark_data::getGenerator().loadReadsFromFastq(fastqFilePath.string(),
                                                             maxReadsToLoad);
}

}  // namespace benchmark_util
