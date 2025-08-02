#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

class Read;     // Forward declare Read
class Options;  // Forward declare options

class Evaluator final {
public:
    explicit Evaluator(Options* opt) noexcept;
    ~Evaluator() = default;

    // Disable copy and enable move, Evaluator is cheap
    Evaluator(const Evaluator&)                    = delete;
    auto operator=(const Evaluator&) -> Evaluator& = delete;
    Evaluator(Evaluator&&)                         = default;
    auto operator=(Evaluator&&) -> Evaluator&      = default;

    auto isTwoColorSystem() -> bool;
    void evaluateSeqLen();
    void evaluateOverRepSeqs();
    void evaluateReadNum(long& readNum);  // evaluate how many reads are stored in the input file
    auto evalAdapterAndReadNum(long& readNum, bool isR2) -> std::string;

    static auto computeSeqLen(const std::string& filename) -> int;
    static void computeOverRepSeq(const std::string&                     filename,
                                  std::unordered_map<std::string, long>& hotseqs,
                                  int                                    seqLen);
    static auto test() -> bool;
    static auto matchKnownAdapter(const std::string& seq) -> std::string;

private:
    auto int2seq(unsigned int val, int seqlen) const -> std::string;
    auto seq2int(const std::string& seq, int pos, int keylen, int lastVal = -1) const -> int;
    auto getAdapterWithSeed(int                                       seed,
                            const std::vector<std::unique_ptr<Read>>& loadedReads,
                            int                                       keylen) -> std::string;
    auto checkKnownAdapters(const std::vector<std::unique_ptr<Read>>& reads) -> std::string;

    Options* mOptions {nullptr};
};
