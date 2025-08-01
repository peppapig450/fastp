#pragma once

#include <memory>
#include <string>
#include <vector>

#include "../../src/read.h"

inline std::vector<std::unique_ptr<Read>> createReadsWithKnownAdapter() {
    const std::string                  adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    std::vector<std::unique_ptr<Read>> reads;
    reads.emplace_back(std::make_unique<Read>("@r1", adapter + "AAAA"));
    reads.emplace_back(std::make_unique<Read>("@r2", std::string("TTTT") + adapter));
    reads.emplace_back(std::make_unique<Read>("@r3", adapter + "CCCC"));
    return reads;
}

inline std::vector<std::unique_ptr<Read>> createReadsWithoutAdapter() {
    std::vector<std::unique_ptr<Read>> reads;
    reads.emplace_back(std::make_unique<Read>("@r1", "ACGTACGTACGT"));
    reads.emplace_back(std::make_unique<Read>("@r2", "TGCATGCATGCA"));
    reads.emplace_back(std::make_unique<Read>("@r3", "GGGGTTTTCCCC"));
    return reads;
}