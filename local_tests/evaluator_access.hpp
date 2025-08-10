#pragma once

#include "evaluator.h"

#include <memory>
#include <string>
#include <vector>

struct EvaluatorAccess {
    static auto checkKnownAdapters(const Evaluator& evaluator,
                                   const std::vector<std::unique_ptr<Read>>& reads) -> std::string {
        return evaluator.checkKnownAdapters(reads);
    }

    static auto seq2int(const Evaluator& evaluator,
                        const std::string& seq,
                        int pos,
                        int keylen,
                        int lastVal = -1) -> int {
        return evaluator.seq2int(seq, pos, keylen, lastVal);
    }

    static auto int2seq(const Evaluator& evaluator, unsigned int val, int seqlen) -> std::string {
        return evaluator.int2seq(val, seqlen);
    }

    static auto getAdapterWithSeed(const Evaluator& evaluator,
                                   int seed,
                                   const std::vector<std::unique_ptr<Read>>& reads,
                                   int keylen) -> std::string {
        return evaluator.getAdapterWithSeed(seed, reads, keylen);
    }
};