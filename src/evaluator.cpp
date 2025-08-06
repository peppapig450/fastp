#include "evaluator.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "fastqreader.h"
#include "knownadapters.h"
#include "nucleotidetree.h"
#include "options.h"
#include "read.h"
#include "util.h"

#ifdef ENABLE_COVERAGE
    #define EXCLUDE_FROM_COVERAGE __attribute__((no_profile_instrument_function))
#else
    #define EXCLUDE_FROM_COVERAGE
#endif

namespace { // anon
namespace detail { // detail

// NEXTSEQ500, NEXTSEQ 550/550DX, NOVASEQ
constexpr std::array<const char*, 4> kTwoColorPrefixes = {"@NS", "@NB", "@NDX", "@A0"};

auto starts_with_any(const std::string& name) -> bool {
    return std::any_of(kTwoColorPrefixes.begin(), kTwoColorPrefixes.end(), [&](const char* prefix) {
        const std::size_t len = strlen(prefix);
        return name.size() >= len && name.compare(0, len, prefix) == 0;
    });
}

// Heuristic for minimum threshold count
auto minCountForLen(int fragmentLen, int seqlen) -> int {
    return (fragmentLen >= seqlen - 1) ? 3
           : (fragmentLen >= 100)      ? 5
           : (fragmentLen >= 40)       ? 20
           : (fragmentLen >= 20)       ? 100
                                       : 500;  // fragmentLen >= 10
}

// Helpers for computeOverRepSeq

using CountMap   = std::unordered_map<std::string, long>;
// Using a pointer allows us to reference the contents of the string without copying
using StringPair = std::pair<const std::string*, long>;

// Count all k-mers at the chosen window sizes for a single read, this is pass 1 for
// computeOverRepSeq
void countKmers(const Read& read, const std::array<int, 5>& steps, CountMap& counts) {
    const auto& seq    = read.seq();
    const auto  seqLen = static_cast<int>(read.length());

    for (int kmerLength : steps) {
        if (seqLen <= kmerLength) {
            continue;  // too short for this k-mer size
        }

        for (int i = 0; i + kmerLength <= seqLen; ++i) {
            ++counts[seq.substr(static_cast<std::size_t>(i), static_cast<std::size_t>(kmerLength))];
        }
    }
}

// Pass 2 for computeOverRepSeq: transfer only entries passing the dynamic threshold into
// `filteredCounts`
void filterByMinCount(const CountMap& inputCounts, CountMap& filteredCounts, int fullSeqLen) {
    using namespace detail;

    for (const auto& countPair : inputCounts) {
        if (countPair.second
            >= minCountForLen(static_cast<int>(countPair.first.length()), fullSeqLen)) {
            filteredCounts.emplace(countPair);
        }
    }
}

// Pass 3 for computeOverRepSeq: prune any short sequence that is largely contained in an longer,
// more abundant one
void eraseContainedSeqs(CountMap& seqs) {
    // Sort pointers by descending length to visit long -> short
    std::vector<StringPair> sorted;
    sorted.reserve(seqs.size());
    for (auto& countPair : seqs) {
        sorted.emplace_back(&countPair.first, countPair.second);
    }

    // Since C++11 does not allow us to use auto in lambda parameters, we use a type alias to avoid
    // typing out the type signature everytime
    std::sort(sorted.begin(), sorted.end(), [](const StringPair& a, const StringPair& b) {
        return a.first->length() > b.first->length();
    });

    const auto sortedCount = sorted.size();

    std::unordered_set<std::string> doomed;
    for (std::size_t i = 0; i < sortedCount; ++i) {
        const auto&        longerPair  = sorted[i];
        const std::string& longer      = *longerPair.first;
        long               longerCount = longerPair.second;

        for (std::size_t j = i + 1; j < sortedCount; ++j) {
            const auto&        shorterPair  = sorted[j];
            const std::string& shorter      = *shorterPair.first;
            long               shorterCount = shorterPair.second;

            if (longer.find(shorter) != std::string::npos && shorterCount < longerCount * 10) {
                doomed.insert(shorter);
            }
        }
    }

    for (const auto& seq : doomed) {
        seqs.erase(seq);
    }
}

// Helpers for checkKnownAdapters

constexpr std::size_t kMaxCheckReads  = 100000;                 // Stop checking after 100,000 reads
constexpr std::size_t kMaxCheckBases  = kMaxCheckReads * 1000;  // Or after 100 million bases
constexpr int         kMaxHit         = 1000;  // Max adapter hit count before saturation
constexpr int         kMinMatch       = 8;     // Minimum required consecutive matching bases
constexpr int         kMismatchFactor = 16;    // Allow 1 mismatch per 16 matched bases

// Helper struct for adapter statistics
struct AdapterStats {
    int hits       = 0;
    int mismatches = 0;
};

using AdapterStatsMap = std::unordered_map<std::string, AdapterStats>;

// Initialize adapter statistics map
auto initAdapterStats(const std::unordered_map<std::string, std::string>& adapters)
    -> AdapterStatsMap {
    AdapterStatsMap stats;
    stats.reserve(adapters.size());

    for (const auto& adapterPair : adapters) {
        stats.emplace(adapterPair.first, AdapterStats {});
    }

    return stats;
}

// Compare adapter to read at a given position
// Returns true if mismatches <= allowed, and sets `mismatches`
auto matchAtPosition(const char* read,
                     int         readLen,
                     const char* adapter,
                     int         adapterLen,
                     int         startPos,
                     int&        mismatches) noexcept -> bool {
    const int compLen         = std::min(readLen - startPos, adapterLen);
    const int allowedMismatch = compLen / kMismatchFactor;
    mismatches                = 0;

    const char* readStart = read + startPos;

    // Loop until we encounter too many mismatches
    for (int i = 0; i < compLen; ++i) {
        if (adapter[i] != readStart[i] && ++mismatches > allowedMismatch) {
            return false;
        }
    }

    return true;
}

// Process a single adapter against a read sequence
void processAdapter(const std::string& seq,
                    const std::string& adapter,
                    AdapterStats&      stats,
                    int&               bestHitSoFar) {
    const int readLen    = static_cast<int>(seq.length());
    const int adapterLen = static_cast<int>(adapter.length());

    // Early exit if it's impossible to fit
    if (adapterLen >= readLen) {
        return;
    }

    // Heuristic to skip weak contestants
    if (bestHitSoFar > 20 && stats.hits < bestHitSoFar / 10) {
        return;
    }

    const char* readData    = seq.data();
    const char* adapterData = adapter.data();
    const int   scanLimit   = readLen - kMinMatch;

    for (int pos = 0; pos < scanLimit; ++pos) {
        int mismatches;
        if (!matchAtPosition(readData, readLen, adapterData, adapterLen, pos, mismatches)) {
            continue;  // no good at this position
        }

        ++stats.hits;
        stats.mismatches += mismatches;
        bestHitSoFar      = std::max(bestHitSoFar, stats.hits);
        break;  // found a hit -> next adapter
    }
}

// Helper to decide which adapter "wins" when we are forced to choose from multiple
auto selectBestAdapter(const AdapterStatsMap& stats, AdapterStats& bestStats)
    -> const std::string* {
    const std::string* bestAdapter = nullptr;

    for (const auto& adapterEntry : stats) {
        const auto& adapterName  = adapterEntry.first;
        const auto& adapterStats = adapterEntry.second;

        const int currentHits       = adapterStats.hits;
        const int currentMismatches = adapterStats.mismatches;
        const int bestHits          = bestStats.hits;
        const int bestMismatches    = bestStats.mismatches;

        if (bestAdapter == nullptr || currentHits > bestHits
            || (currentHits == bestHits && currentMismatches < bestMismatches)
            || (currentHits == bestHits && currentMismatches == bestMismatches
                && adapterName.length() < bestAdapter->length())) {
            bestAdapter = &adapterName;
            bestStats   = adapterStats;
        }
    }

    return bestAdapter;
}

// Simple confidence heuristic helper
auto isConfidentMatch(const AdapterStats& stats, std::size_t checkedReads) -> bool {
    return (stats.hits > static_cast<int>(checkedReads / 50))
           || (stats.hits > static_cast<int>(checkedReads / 200)
               && stats.mismatches < static_cast<int>(checkedReads));
}

// Helpers for evalAdapterAndReadNum

// Minimum reads required for evaluation
constexpr int kMinRecordsEval    = 10000;
// Length of k-mer
constexpr int kKeyLen            = 10;
// Total number of unique k-mers (4^k)
constexpr std::size_t kKeySpace  = 1ULL << (kKeyLen * 2ULL);
// Number of top scoring k-mers to retain
constexpr int kTopCandidates     = 10;
// Minimum fold-enrichment to consider a k-mer significant
constexpr int kFoldThreshold     = 20;
// We use a 1% safety limit to account for under-evaluation due to potential bad quality
constexpr double kSafetyFactor = 1.01;

struct SampleResult {
public:
    std::vector<std::unique_ptr<Read>> reads;  // sampled reads (owns memory if collectSeq=true)
    std::size_t                        readsCount  = 0;      // number of reads sampled
    std::size_t                        bases       = 0;      // total bases sampled
    std::size_t                        bytesRead   = 0;      // current bytes position in file
    std::size_t                        bytesTotal  = 0;      // total bytes in file
    std::size_t                        firstOffset = 0;      // byte offset of 1st record
    bool                               reachedEOF  = false;  // true -> consumed whole file

    // Linear extrapolation : bytes / reads => total reads
    auto extrapolateTotalReads(double safety = kSafetyFactor) const -> long {
        if (readsCount == 0) {
            return 0L;
        }

        if (reachedEOF) {
            return static_cast<long>(readsCount);
        }

        const double bytesPerRead =
            static_cast<double>(bytesRead - firstOffset) / static_cast<double>(readsCount);

        return static_cast<long>(static_cast<double>(bytesTotal) * safety / bytesPerRead);
    }
};

/**
 * Core sampler that powers all sampling needs in `Evaluator`.
 *
 * @param fastqPath    Input FASTQ.
 * @param maxReads     Stop sampling after this many reads.
 * @param maxBases     Stop sampling after this many bases.
 * @param collectSeq   When *true* the sequences are kept in memory (needed by
 *                     de-novo adapter discovery); otherwise only statistics
 *                     are gathered.
 */
auto sampleFastq(const std::string& fastqPath,
                 std::size_t        maxReads,
                 std::size_t        maxBases,
                 bool               collectSeq = false) -> SampleResult {
    SampleResult result;
    FastqReader  reader {fastqPath};

    while (result.readsCount < maxReads && result.bases < maxBases) {
        std::unique_ptr<Read> read {reader.read()};
        if (read == nullptr) {
            result.reachedEOF = true;
            break;
        }

        // We sample the file pointer for byte-offsets on the first read
        if (result.readsCount == 0) {
            reader.getBytes(result.bytesRead, result.bytesTotal);
            result.firstOffset = result.bytesRead;
        }

        const std::size_t readLen  = read->length();
        result.bases              += readLen;
        ++result.readsCount;

        if (collectSeq) {
            result.reads.emplace_back(std::move(read));
        }
    }

    // Refresh the byte counter if we don't reach the end of the file
    if (!result.reachedEOF && result.readsCount > 0) {
        reader.getBytes(result.bytesRead, result.bytesTotal);
    }

    return result;
}

struct SeedCandidate {
    int          key   = -1;  // packed 2-bit representation of the k-mer
    unsigned int count = 0;   // how often it was observed

    SeedCandidate() = default;
    explicit SeedCandidate(int key_, unsigned int count_)
        : key(key_)
        , count(count_) {}
};

// Quick low-complexity screen
auto isLowComplexityKey(int packedKey) -> bool {
    // clang-format off
    std::array<int, 4> baseCount{{0,0,0,0,}};
    // clang-format on
    auto key = static_cast<unsigned int>(packedKey);

    for (int i = 0; i < kKeyLen; ++i) {
        ++baseCount[(key >> (i*2U)) & 0x3U];
    }

    if (*std::max_element(baseCount.begin(), baseCount.end()) >= kKeyLen - 4) {
        return true;
    }

    return (baseCount[2] + baseCount[3] >= kKeyLen - 2) || ((key >> 12U) == 0xFFU);
}

// Pick the top-N most enriched seeds that survive complexity / enrichment filtering
auto selectTopSeeds(const std::vector<unsigned int>& kmerCounts, std::uint64_t totalCounts)
    -> std::vector<SeedCandidate> {
    std::vector<SeedCandidate> topSeeds;
    topSeeds.reserve(kTopCandidates);

    const std::uint64_t enrichmentCutoff =
        (totalCounts * static_cast<std::uint64_t>(kFoldThreshold)) / kKeySpace;

    for (std::size_t k = 0; k < kmerCounts.size(); ++k) {
        const unsigned int count = kmerCounts[k];
        if (count < 10) {
            continue;
        }

        if (isLowComplexityKey(static_cast<int>(k))) {
            continue;
        }

        if (count < enrichmentCutoff) {
            continue;
        }

        // insert into descending ordered container
        SeedCandidate candidate {static_cast<int>(k), count};

        auto iter = std::upper_bound(topSeeds.begin(),
                                     topSeeds.end(),
                                     candidate,
                                     [](const SeedCandidate& a, const SeedCandidate& b) -> bool {
                                         return a.count > b.count;
                                     });
        topSeeds.insert(iter, candidate);
        if (topSeeds.size() > kTopCandidates) {
            topSeeds.pop_back();
        }
    }

    return topSeeds;
}

}  // namespace detail
}  // namespace

Evaluator::Evaluator(Options* opt) noexcept : mOptions(opt) {}

bool Evaluator::isTwoColorSystem() {
    using namespace detail;

    FastqReader reader(mOptions->in1);

    // Wrap in a unique_ptr for now
    std::unique_ptr<Read> r{reader.read()}; // TODO: update fastqreader to return a unique ptr

    if (r == nullptr) {
        return false;
    }

    // NEXTSEQ500, NEXTSEQ 550/550DX, NOVASEQ
    const std::string& name = r->name();
    return starts_with_any(name);
}

void Evaluator::evaluateSeqLen() {
    const auto& in1 = mOptions->in1;
    const auto& in2 = mOptions->in2;

    if (!in1.empty()) {
        mOptions->seqLen1 = computeSeqLen(in1);
    }
    if (!in2.empty()) {
        mOptions->seqLen2 = computeSeqLen(in2);
    }
}

int Evaluator::computeSeqLen(const std::string& filename) {
    FastqReader           reader(filename);
    constexpr std::size_t kRecordLimit = 1000;
    int                   maxLen       = 0;

    for (std::size_t rec = 0; rec < kRecordLimit; ++rec) {
        std::unique_ptr<Read> r {reader.read()};
        if (r == nullptr) {
            break;
        }

        maxLen = std::max(maxLen, static_cast<int>(r->length()));
    }

    return maxLen;
}

void Evaluator::computeOverRepSeq(const std::string&                filename,
                                  std::unordered_map<string, long>& hotseqs,
                                  int                               seqlen) {
    using namespace detail;

    constexpr std::int64_t   kBaseLimit = 151LL * 10000LL;  // 1.51 M bases ≈ 10k 150 bp reads
    const std::array<int, 5> kSteps {{10, 20, 40, 100, std::min(150, seqlen - 2)}};

    FastqReader reader(filename);

    // We reserve a decent-sized bucket count to reduce re-hashing during counting
    CountMap allCounts;
    // TODO: this is a heuristic and should be tuned to sample size, since sample sizes can vary we
    // should probably do this at runtime, maybe using the getBytes method
    allCounts.reserve(500000);

    std::int64_t basesSeen = 0;

    // First pass we count k-mers
    while (basesSeen < kBaseLimit) {
        std::unique_ptr<Read> rec(reader.read());
        if (rec == nullptr) {
            break;
        }

        const Read& read  = *rec;
        basesSeen        += static_cast<std::int64_t>(read.length());
        countKmers(read, kSteps, allCounts);
    }

    // Move frequent fragments into caller-supplied map then prune
    filterByMinCount(allCounts, hotseqs, seqlen);
    eraseContainedSeqs(hotseqs);

    // output for test
    /*for(iter = hotseqs.begin(); iter!=hotseqs.end(); iter++) {
        cerr << iter->first << ": " << iter->second << endl;
    }*/
}

void Evaluator::evaluateOverRepSeqs() {
    auto& opts = *mOptions;

    const auto& in1 = opts.in1;
    const auto& in2 = opts.in2;

    if (!in1.empty()) {
        computeOverRepSeq(in1, opts.overRepSeqs1, opts.seqLen1);
    }

    if (!in2.empty()) {
        computeOverRepSeq(in2, opts.overRepSeqs2, opts.seqLen2);
    }
}

void Evaluator::evaluateReadNum(long& readNum) {
    using namespace detail;

    // Limit reads to ~0.5 million records
    constexpr std::int64_t kReadLimit = 512LL * 1024;
    // Same cap in bases (NOTE: This assumes 150bp that may not be optimal,
    // but for consistency with the original code we keep it)
    constexpr std::int64_t kBaseLimit = 151LL * kReadLimit;

    const auto stats = sampleFastq(mOptions->in1, kReadLimit, kBaseLimit, /*collectSeq=*/false);
    readNum          = stats.extrapolateTotalReads();
}

auto Evaluator::checkKnownAdapters(const ReadsVector& reads) const -> std::string {
    using namespace detail;

    const auto& knownAdapters = adapters::getKnown();
    auto        stats         = initAdapterStats(knownAdapters);

    std::size_t checkedReads = 0;
    std::size_t checkedBases = 0;
    int         bestHitSoFar = 0;

    for (const auto& readPtr : reads) {
        const auto& seq     = readPtr->seq();
        const int   readLen = static_cast<int>(seq.length());

        ++checkedReads;
        checkedBases += static_cast<std::size_t>(readLen);

        if (checkedReads > kMaxCheckReads || checkedBases > kMaxCheckBases
            || bestHitSoFar > kMaxHit) {
            break;  // Exit once we have enough evidence
        }

        // Compare against every known adapter
        for (const auto& adapterPair : knownAdapters) {
            auto it = stats.find(adapterPair.first);
            processAdapter(seq, adapterPair.first, it->second, bestHitSoFar);
        }
    }

    // Choose the best fit adapter
    AdapterStats       winnerStats;
    const std::string* bestAdapter = selectBestAdapter(stats, winnerStats);
    if (bestAdapter == nullptr || !isConfidentMatch(winnerStats, checkedReads)) {
        return {};  // No clear adapter found
    }

    auto chosenAdapter = *bestAdapter;
    std::cerr << knownAdapters.at(chosenAdapter) << '\n' << chosenAdapter << '\n';
    return chosenAdapter;
}

auto Evaluator::evalAdapterAndReadNum(long& readNum, bool isR2) -> std::string {
    using namespace detail;

    // Max number of reads to sample (<= 256,000)
    constexpr std::size_t kReadLimit = 256ULL * 1024ULL;
    // Upper limit of bases assuming 151 bp reads
    // NOTE: since sequencers such as Illumina MiSeq i100 support 300 bp reads this may be too small
    constexpr std::size_t kBaseLimit = 151ULL * kReadLimit;

    const std::string& fastqPath = isR2 ? mOptions->in2 : mOptions->in1;

    // sample and extrapolate read count
    auto reads = sampleReadsForAdapter(fastqPath, readNum, kReadLimit, kBaseLimit);

    // try known adapter list first
    std::string known = checkKnownAdapters(reads);
    if (!known.empty() && known.length() > static_cast<std::size_t>(kMinMatch)) {
        return known;
    }

    // final fallback to de-novo inference (seed enrichment + extension)
    return inferAdapterDeNovo(reads);
}

auto Evaluator::sampleReadsForAdapter(const std::string& fastqPath,
                                      long&              readNum,
                                      std::size_t        maxReads,
                                      std::size_t        maxBases) const -> ReadsVector {
    using namespace detail;

    // TODO: we may want to explicitly std::move here since RVO is not guaranteed
    auto sample = sampleFastq(fastqPath, maxReads, maxBases, /*collectSeq=*/true);

    readNum = sample.extrapolateTotalReads();

    // unique_ptr is not copyable so this MUST be moved
    return std::move(sample.reads);
}

auto Evaluator::hasSufficientSample(const ReadsVector& reads) const -> bool {
    using namespace detail;
    return reads.size() >= static_cast<std::size_t>(kMinRecordsEval);
}

auto Evaluator::inferAdapterDeNovo(const ReadsVector& reads) const -> std::string {
    using namespace detail;

    // build k-mer histogram, skipping noisy tail data common in Illumina data with shift)
    const std::size_t         shiftTail = Evaluator::tailShift();
    std::vector<unsigned int> kmerCounts(kKeySpace, 0);
    const std::uint64_t       totalCounts = countKmersForAdapter(reads, kmerCounts, shiftTail);

    // pick enriched, non-low-complexity seeds
    auto seeds = selectTopSeeds(kmerCounts, totalCounts);

    // extend each seed into a candidate adapter, returning the first confident one
    for (const auto& seed : seeds) {
        std::string adapter = getAdapterWithSeed(seed.key, reads, kKeyLen);
        if (!adapter.empty()) {
            return adapter;
        }
    }

    return {};  // No adapter found
}

auto Evaluator::countKmersForAdapter(const ReadsVector&         reads,
                                     std::vector<unsigned int>& kmerCounts,
                                     std::size_t                shiftTail) const-> std::uint64_t {
    using namespace detail;

    std::uint64_t totalCounts = 0;

    constexpr std::size_t kStartPos = 20;

    for (const auto& read : reads) {
        const std::string& seq = read->seq();
        // Skip if the sequence is too short
        if (seq.length() < kStartPos + shiftTail + static_cast<std::size_t>(kKeyLen)) {
            continue;
        }

        const std::size_t limit   = seq.length() - shiftTail - static_cast<std::size_t>(kKeyLen);
        int               rolling = -1;

        for (std::size_t pos = kStartPos; pos <= limit; ++pos) {
            rolling = Evaluator::seq2int(seq, static_cast<int>(pos), kKeyLen, rolling);
            if (rolling >= 0) {
                ++kmerCounts[static_cast<std::size_t>(rolling)];
                ++totalCounts;
            }
        }
    }

    kmerCounts[0] = 0;  // ignore poly-A seed "AAAAAAAAAA"
    return totalCounts;
}

auto Evaluator::getAdapterWithSeed(int seed, const ReadsVector& reads, int keylen) const
    -> std::string {
    constexpr int kMinScanPos       = 20;   // skip low-quality head
    constexpr int kMaxSearchLength  = 500;  // avoid wasting time deep in the read
    constexpr int kMaxAdapterLength = 60;   // hard-cap on reported adapter
    // We shift the last cycle for evaluation since Illumina data tends to be noisy
    const int shiftTail             = std::max(1, mOptions->trim.tail1);

    // Sequence after the seed
    NucleotideTree forwardTree {mOptions};
    // Sequence before the seed (reversed)
    NucleotideTree backwardTree {mOptions};

    const int keyShiftTotal = keylen + shiftTail;
    for (const auto& rec : reads) {
        const auto readLen = static_cast<int>(rec->length());
        const int  limit   = readLen - keyShiftTotal;

        // Skip if the read is too short to be useful
        if (limit <= kMinScanPos) {
            continue;
        }

        const std::string& seq        = rec->seq();
        int                rollingKey = -1;  // reset rolling hash for this read

        for (int pos = kMinScanPos; pos <= limit && pos < kMaxSearchLength; ++pos) {
            rollingKey = seq2int(seq, pos, keylen, rollingKey);
            if (rollingKey != seed) {
                continue;
            }

            const auto spos    = static_cast<std::size_t>(pos);
            const auto skeylen = static_cast<std::size_t>(keylen);
            const auto slimit  = static_cast<std::size_t>(limit);

            // Forward context (suffix, same orientation)
            forwardTree.addSeq(seq.substr(spos + skeylen, slimit - spos));

            // backward context (prefix, reverse orientation)
            backwardTree.addSeq(reverse(seq.substr(0, spos)));
        }
    }

    // Build candidate adapter (5'-prefix + seed + 3'-suffix)
    bool reachedLeafForward  = true;
    bool reachedLeafBackward = true;

    const std::string forwardPath  = forwardTree.getDominantPath(reachedLeafForward);
    const std::string backwardPath = backwardTree.getDominantPath(reachedLeafBackward);

    std::string adapter = reverse(backwardPath)    // re-orient prefix
                          + int2seq(seed, keylen)  // seed itself
                          + forwardPath;           // suffix

    if (adapter.length() > kMaxSearchLength) {
        adapter.resize(kMaxAdapterLength);
    }

    // Prefer a perfect match against the adapter list
    auto matched = adapters::matchKnown(adapter);
    if (!matched.empty()) {
        const auto& known = adapters::getKnown();
        std::cerr << known.at(matched) << '\n' << matched << '\n';
        return matched;
    }

    // Otherwise return the de-novo adapter only if both paths were resolved
    if (reachedLeafForward && reachedLeafBackward) {
        std::cerr << adapter << '\n';
        return adapter;
    }

    return {};
}

auto Evaluator::int2seq(std::uint32_t val, int seqlen) const -> std::string {
    // Mapping for 2-bit values to DNA bases:
    // 00 -> 'A', 01 -> 'T', 10 -> 'C', 11 -> 'G'
    static constexpr std::array<char, 4> kBases {{'A', 'T', 'C', 'G'}};

    // Create a string with length `seqlen`, initialized with 'N' (unknown base)
    std::string result(seqlen, 'N');

    // Convert the integer to a sequence of DNA bases (2 bits per base)
    for (int i = seqlen - 1; i >= 0; --i) {
        // Extract the lowest 2 bits of `val` (mask with 0b11 = 0x3)
        // These 2 bits represent one DNA base (since 2 bits can encode 4 values)
        result[i] =
            kBases[val & 0x3U];  // TODO: if 0x3U is used elsewhere make it a file level constant

        // Shift `val` right by 2 bases to process the next base in the next iteration
        val >>= 2;
    }

    return result;
}

/*
 * Converts a DNA substing (A/T/C/G) to a packed 2-bit integer representation
 *
 * If `lastVal` >= 0, performs a rolling hash by reusing the previous k-mer value
 * and shifting in one new base. Otherwise, computes the hash from scratch.
 *
 * Each base is encoded as:
 *   A/a -> 00, T/t -> 01, C/c -> 10, G/g -> 11
 *
 * Returns -1 if any character is not a valid base (e.g., 'N').
 */
auto Evaluator::seq2int(const std::string& seq, int pos, int keylen, int lastVal) const -> int {
    constexpr int kInvalid = -1;

    // Map ASCII characters to their 2-bit base encoding, this lookup table keeps the inner loop
    // branch-free and case-insensitive
    static const std::array<int, 256> kBaseToBits = [kInvalid] {
        std::array<int, 256> table {};

        table.fill(kInvalid);
        // clang-format off
        table['A'] = table['a'] = 0;
        table['T'] = table['t'] = 1;
        table['C'] = table['c'] = 2;
        table['G'] = table['g'] = 3;
        // clang-format on

        return table;
    }();

    const auto decode_base = [&](char base) noexcept -> int {
        return kBaseToBits[static_cast<unsigned char>(base)];
    };

    // Perform a rolling update reusing the previous value
    if (lastVal >= 0) {
        // Create a mask to retain only the lowest 2 * keylen bits
        const std::uint32_t bitmask = (1U << (static_cast<std::uint32_t>(keylen) * 2U)) - 1U;

        // Shift previous value left by 2 bits, which drops the oldest base and apply the mask
        std::uint32_t newKey = (static_cast<std::uint32_t>(lastVal) << 2U) & bitmask;

        // Decode any new incoming base at the end of the k-mer window
        const auto baseBits =
            decode_base(seq[pos + keylen - 1]);  // TODO: These should probably be std::size_t
        if (baseBits == kInvalid) {
            return kInvalid;
        }

        // Add the new base to the right end of the key
        return static_cast<int>(newKey | static_cast<std::uint32_t>(baseBits));
    }

    // Compute a new fresh hash from scratch
    std::uint32_t key = 0;
    for (int i = 0; i < keylen; ++i) {
        const auto baseBits = decode_base(seq[pos + i]);
        if (baseBits == kInvalid) {
            return kInvalid;
        }

        // Shift current key by 2 bits and add new base
        key = (key << 2U) | static_cast<std::uint32_t>(baseBits);
    }

    return static_cast<int>(key);
}

// We shift cycles at the end to avoid noisy tail data that is common in Illumina
auto Evaluator::tailShift() const -> std::size_t {
    assert(mOptions->trim.tail1 >= 0);
    return std::max<std::size_t>(std::size_t(1), static_cast<std::size_t>(mOptions->trim.tail1));
}

EXCLUDE_FROM_COVERAGE
bool Evaluator::test() {
    Evaluator eval(nullptr);
    bool passedTests = true;

    // round-trip several sequences through seq2int/int2seq
    std::vector<std::string> seqs = {"ATCGATCGAT", "GGGGGGGGGG", "TATATATATA"};
    for (const auto& s : seqs) {
        int val = eval.seq2int(s, 0, static_cast<int>(s.length()), -1);
        if (eval.int2seq(val, static_cast<int>(s.length())) != s) {
            std::cerr << "round-trip failed for " << s << "\n";
            passedTests = false;
        }
    }

    // verify rolling seq2int produces the same result as computing from scratch
    std::string rollingSeq = "ATCGATCG";
    int keylen = 4;
    int rolling = -1;
    for (int i = 0; i <= static_cast<int>(rollingSeq.length()) - keylen; ++i) {
        rolling = eval.seq2int(rollingSeq, i, keylen, rolling);
        int fromScratch = eval.seq2int(rollingSeq, i, keylen, -1);
        if (rolling != fromScratch) {
            std::cerr << "rolling seq2int mismatch at pos " << i << "\n";
            passedTests = false;
        }
    }

    // seq2int should flag sequences containing invalid bases
    if (eval.seq2int("ATCN", 0, 4, -1) != -1) {
        std::cerr << "seq2int should return -1 for sequences containing N" << "\n";
        passedTests = false;
    }

    // simple int2seq checks
    if (eval.int2seq(0, 3) != "AAA") {
        std::cerr << "int2seq failed for value 0" << "\n";
        passedTests = false;
    }
    if (eval.int2seq((1 << 6) - 1, 3) != "GGG") {
        std::cerr << "int2seq failed for max value" << "\n";
        passedTests = false;
    }

    return passedTests;
}
