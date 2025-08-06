#pragma once

#include <cstddef>
#include <cstdint>
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

private:
    using ReadsVector = std::vector<std::unique_ptr<Read>>;

    auto int2seq(unsigned int val, int seqlen) const -> std::string;
    auto seq2int(const std::string& seq, int pos, int keylen, int lastVal = -1) const -> int;
    auto getAdapterWithSeed(int seed, const ReadsVector& reads, int keylen) const -> std::string;
    auto checkKnownAdapters(const ReadsVector& reads) const -> std::string;

    auto countKmersForAdapter(const ReadsVector&         reads,
                              std::vector<unsigned int>& kmerCounts,
                              std::size_t                shiftTail) const -> std::uint64_t;

    auto tailShift() const -> std::size_t;

    // Helpers for evalAdapterAndReadNum

    // TODO: determine whether these should go into anon namespace or not, the main reason I am not
    // sure if they should is that it would require duplicating ReadsVector, and the functions arent
    // exactly very generic
    auto sampleReadsForAdapter(const std::string& fastqPath,
                               long&              readNum,
                               std::size_t        maxReads,
                               std::size_t        maxBases) const -> ReadsVector;

    auto hasSufficientSample(const ReadsVector& reads) const -> bool;

    auto inferAdapterDeNovo(const ReadsVector& reads) const -> std::string;

    Options* mOptions {nullptr};
};
