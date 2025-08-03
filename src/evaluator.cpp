#include "evaluator.h"

#include <algorithm>
#include <cstdint>
#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
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

// NEXTSEQ500, NEXTSEQ 550/550DX, NOVASEQ
constexpr std::array<const char*, 4> prefixes = {"@NS", "@NB", "@NDX", "@A0"};

auto starts_with_any(const std::string& name) -> bool {
    return std::any_of(prefixes.begin(), prefixes.end(), [&](const char* prefix) {
        std::size_t len = strlen(prefix);
        return name.size() >= len && name.compare(0, len, prefix) == 0;
    });
}
}  // namespace

Evaluator::Evaluator(Options* opt) noexcept : mOptions(opt) {}

bool Evaluator::isTwoColorSystem() {
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

// TODO: This should be split up into multiple functions as it is essentially 3 independent operations
void Evaluator::computeOverRepSeq(const std::string&                filename,
                                  std::unordered_map<string, long>& hotseqs,
                                  int                               seqlen) {
    constexpr std::int64_t   kBaseLimit = 151LL * 10000LL;  // 1.51 M bases ≈ 10 k 150 bp reads
    const std::array<int, 5> kSteps {{10, 20, 40, 100, std::min(150, seqlen - 2)}};

    FastqReader reader(filename);

    // We reserve a decent-sized bucket count to reduce re-hashing during counting
    std::unordered_map<std::string, long> seqCounts;
    // TODO: this is a heuristic and should be tuned to sample size, since sample sizes can vary we
    // should probably do this at runtime, maybe using the getBytes method
    seqCounts.reserve(500000);

    int basesSeen = 0;

    // First pass we count k-mers
    while (basesSeen < kBaseLimit) {
        std::unique_ptr<Read> read(reader.read());
        if (read == nullptr) {
            break;
        }

        const auto& seq     = read->seq();
        const int   readLen = static_cast<int>(seq.size());
        basesSeen += readLen;

        for (int step : kSteps) {
            if (readLen <= step) {
                continue;  // too short for this k
            }

            for (int i = 0; i + step <= readLen; ++i) {
                // operator[] avoids double lookup
                seqCounts[seq.substr(i, step)]++;
            }
        }
    }

    // Apply length based thresholds
    const auto minCount = [seqlen](int len) -> int {
        return (len >= seqlen - 1) ? 3
               : (len >= 100)      ? 5
               : (len >= 40)       ? 20
               : (len >= 20)       ? 100
                                   : 500;  // len >= 10
    };

    // because unordered_map uses const's internally, and extract doesn't exist yet (C++17),
    // this is going to always copy the string
    for (const auto& kv : seqCounts) {
        if (kv.second >= minCount(static_cast<int>(kv.first.size()))) {
            hotseqs.emplace(kv);
        }
    }

    // Since C++11 does not allow us to use auto in lambda parameters, we use a type alias to avoid
    // typing out the type signature everytime
    using PtrCountPair = std::pair<const std::string*, long>;

    // we use pointers here to string so that we avoid copying string objects directly,
    // we can manipulate/move the original string.
    std::vector<PtrCountPair> sorted;
    sorted.reserve(hotseqs.size());
    for (auto& kv : hotseqs) {
        // Take a pointer to the string key
        sorted.emplace_back(&kv.first, kv.second);
    }

    std::sort(sorted.begin(), sorted.end(), [](const PtrCountPair& a, const PtrCountPair& b) {
        return a.first->size() > b.first->size();
    });

    // remove substrings
    std::unordered_set<std::string> toErase;
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        const auto&        longerPair  = sorted[i];
        const std::string& longer      = *longerPair.first;
        long               longerCount = longerPair.second;

        for (std::size_t j = i + 1; j < sorted.size(); ++j) {
            const auto&        shorterPair  = sorted[j];
            const std::string& shorter      = *shorterPair.first;
            long               shorterCount = shorterPair.second;

            if (longer.find(shorter) != std::string::npos && shorterCount < longerCount * 10) {
                toErase.insert(shorter);
            }
        }
    }

    for (const auto& seq : toErase) {
        hotseqs.erase(seq);
    }

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
    // Limit reads to ~0.5 million records
    constexpr std::int64_t kReadLimit = 512LL * 1024;
    // Same cap in bases (NOTE: This assumes 150bp that may not be optimal,
    // but for consistency with the original code we keep it)
    constexpr std::int64_t kBaseLimit = 151LL * kReadLimit;
    // +1% head room
    constexpr double kSafetyFactor    = 1.01;

    FastqReader reader {mOptions->in1};

    std::size_t bytesRead      = 0;
    std::size_t bytesTotal     = 0;
    std::size_t firstRecOffset = 0;

    std::int64_t sampledReads = 0;
    std::int64_t sampledBases = 0;
    bool         reachedEOF   = false;

    while (sampledReads < kReadLimit && sampledBases < kBaseLimit) {
        std::unique_ptr<Read> rec {reader.read()};
        // If the returned pointer is null, there are no reads left so we mark EOF and break
        if (rec == nullptr) {
            reachedEOF = true;
            break;
        }

        // After the very first read we get the entire file size
        if (sampledReads == 0) {
            reader.getBytes(bytesRead, bytesTotal);
            firstRecOffset = bytesRead;
        }

        ++sampledReads;
        sampledBases += rec->length();
    }

    readNum = 0;
    // For small files where we reach the end of the file within the limit we use the exact count
    if (reachedEOF) {
        readNum = sampledReads;
        // For larger files exceeding the read limit we extrapolate progress so far
    } else if (sampledReads > 0) {
        // We call getBytes to update bytesRead so we don't need to re-evaluate if splitting output
        // is enabled
        reader.getBytes(bytesRead, bytesTotal);
        const double bytesPerRead =
            static_cast<double>(bytesRead - firstRecOffset) / static_cast<double>(sampledReads);

        // We use a 1% safety limit to account for under-evaluation due to potential bad quality
        readNum = static_cast<long>(static_cast<double>(bytesTotal) * kSafetyFactor / bytesPerRead);
    }
}

// TODO: This should be split up into multiple functions
std::string Evaluator::checkKnownAdapters(const std::vector<std::unique_ptr<Read>>& reads) {
    const auto& knownAdapters = adapters::getKnown();

    // TODO: giving Stats some methods for repeated code might be good
    struct Stats {
        int hits       = 0;  // number of matching reads
        int mismatches = 0;  // total mismatched bases among hits
    };

    // Per-adapter statistics, (we reserve to avoid rehashing)
    std::unordered_map<std::string, Stats> stats;
    stats.reserve(knownAdapters.size());
    for (const auto& kv : knownAdapters) {
        stats.emplace(kv.first, Stats {});
    }

    constexpr std::size_t kMaxCheckReads  = 100000;                 // Allow up to 100k reads
    constexpr std::size_t kMaxCheckBases  = kMaxCheckReads * 1000;  // Allow up to 100M bases
    constexpr int         kMaxHit         = 1000;                   // at 1000 hits we exit
    constexpr int         kMinMatch       = 8;   // minimal consecutive bases to compare
    constexpr int         kMismatchFactor = 16;  // allow 1 mismatch every N bases

    std::size_t checkedReads = 0;
    std::size_t checkedBases = 0;
    int         bestHitSoFar = 0;  // running best hit count across adapters

    for (const auto& readPtr : reads) {
        const auto& seq     = readPtr->seq();
        const auto  readLen = static_cast<int>(seq.length());

        ++checkedReads;
        checkedBases += static_cast<std::size_t>(readLen);
        if (checkedReads > kMaxCheckReads || checkedBases > kMaxCheckBases
            || bestHitSoFar > kMaxHit) {
            break;  // Exit once we have enough evidence
        }

        const char* const readData = seq.data();

        // Try every known adapter against this read
        for (const auto& adapterPair : knownAdapters) {
            const auto& adapter = adapterPair.first;
            auto&       st      = stats[adapter];

            const auto adapterLen = static_cast<int>(adapter.length());

            // Skip if the adapter is longer than the read, as it's impossible to match
            if (adapterLen >= readLen) {
                continue;
            }

            // Heuristic: skip unlikely adapters to save work
            if (bestHitSoFar > 20 && st.hits < bestHitSoFar / 10) {
                continue;
            }

            const char* const adapterData = adapter.data();
            const int         scanLimit   = readLen - kMinMatch;

            // NOTE: this logic is pretty much identical to that in adaptertrimmer.cpp
            for (int pos = 0; pos < scanLimit; ++pos) {
                const int compLen         = std::min(readLen - pos, adapterLen);
                const int allowedMismatch = compLen / kMismatchFactor;
                int       mismatch        = 0;

                const char* readPtr = &readData[pos];
                for (int i = 0; i < compLen; ++i) {
                    // Abort the position if there are too many mismatches
                    if (adapterData[i] != readPtr[i] && ++mismatch > allowedMismatch) {
                        break;
                    }
                }

                if (mismatch <= allowedMismatch) {
                    ++st.hits;
                    st.mismatches += mismatch;
                    bestHitSoFar   = std::max(bestHitSoFar, st.hits);
                    break;  // Match found, we're done with this adapter
                }
            }
        }
    }

    // Choose the best fit adapter: more hits -> fewer mismatches -> shorter length
    const std::string* bestAdapter = nullptr;
    Stats              bestStats {};

    for (const auto& kv : stats) {
        const auto& name = kv.first;
        const auto& st   = kv.second;

        const int currentHits       = st.hits;
        const int currentMismatches = st.mismatches;
        const int bestHits          = bestStats.hits;
        const int bestMismatches    = bestStats.mismatches;

        if (bestAdapter == nullptr || currentHits > bestHits
            || (currentHits == bestHits && currentMismatches < bestMismatches)
            || (currentHits == bestHits && currentMismatches == bestMismatches
                && name.length() < bestAdapter->length())) {
            bestAdapter = &name;
            bestStats   = st;
        }
    }

    // No evidence of any adapter at all
    if (bestAdapter == nullptr) {
        return {};
    }

    const int  maxHits   = bestStats.hits;
    const bool confident = (maxHits > static_cast<int>(checkedReads / 50))
                           || (maxHits > static_cast<int>(checkedReads / 200)
                               && bestStats.mismatches < static_cast<int>(checkedReads));

    if (confident) {
        const auto chosenAdapter = *bestAdapter;
        std::cerr << knownAdapters.at(chosenAdapter) << '\n' << chosenAdapter << '\n';
        return chosenAdapter;
    }

    return {};
}

// TODO: a lot of this logics is duplicated from elsewhere in the file, and it needs to be split up as well
std::string Evaluator::evalAdapterAndReadNum(long& readNum, bool isR2) {
    constexpr std::size_t kReadLimit      = 256ULL * 1024ULL;  // 256k reads
    // Same cap in bases (~= readLimit * 151bp)
    constexpr std::size_t kBaseLimit      = 151ULL * kReadLimit;
    constexpr int         kKeyLen         = 10;  // k-mer length used for adapter seed
    constexpr std::size_t kKeySpace       = 1ULL << (kKeyLen * 2ULL);  // 4^kKeyLen possible k-mers
    constexpr int         kTopCandidates  = 10;                        // keep top-N enriched k-mers
    constexpr int         kFoldThreshold  = 20;     // enrichment fold for candidate selection
    constexpr int         kMinRecordsEval = 10000;  // need at least 10k reads to continue
    // TODO: this is duplicated, move it to anon namesapce
    constexpr double kSafetyFactor        = 1.01;  // +1% head-room when extrapolating

    const auto& filename = isR2 ? mOptions->in2 : mOptions->in1;
    FastqReader reader {filename};

    std::vector<std::unique_ptr<Read>> reads;
    reads.reserve(kReadLimit);

    std::size_t bytesRead      = 0;
    std::size_t bytesTotal     = 0;
    std::size_t firstRecOffset = 0;

    std::size_t sampledBases = 0;
    bool        reachedEOF   = false;

    while (reads.size() < kReadLimit && sampledBases < kBaseLimit) {
        std::unique_ptr<Read> rec {reader.read()};

        if (rec == nullptr) {
            reachedEOF = true;
            break;
        }

        // mark the start of the first record so that we can estimate bytes/read
        if (reads.empty()) {
            reader.getBytes(bytesRead, bytesTotal);
            firstRecOffset = bytesRead;
        }

        sampledBases += rec->length();
        reads.emplace_back(std::move(rec));
    }

    readNum = 0;
    if (reachedEOF) {
        readNum = static_cast<long>(reads.size());
    } else if (!reads.empty()) {
        // update readNum so we avoid having to re-evaluate if splitting output is enabled
        reader.getBytes(bytesRead, bytesTotal);
        const double bytesPerRead =
            static_cast<double>(bytesRead - firstRecOffset) / static_cast<double>(reads.size());
        // We use a 1% safety limit to account for under-evaluation due to potential bad quality
        readNum = static_cast<long>(static_cast<double>(bytesTotal) * kSafetyFactor / bytesPerRead);
    }

    // Early exit, returning empty string if there is not enough data to evaluate reliably
    if (reads.size() < kMinRecordsEval) {
        return {};
    }

    // Try to match against the known adapter list
    auto known = checkKnownAdapters(reads);
    if (known.length() > 8) {
        return known;
    }

    // Shift last cycle(s) to avoid noisy, low-quality tail bases (esp. in Illumina); uses
    // trim.tail1
    const std::size_t shiftTail = std::max<std::size_t>(1, mOptions->trim.tail1);

    std::vector<unsigned int> kmerCounts(kKeySpace, 0);

    for (const auto& rec : reads) {
        const std::string& seq = rec->seq();
        if (seq.length() < 20 + shiftTail + kKeyLen) {
            continue;  // read too short for k-mer scan
        }

        const std::size_t limit      = seq.length() - shiftTail - kKeyLen;
        int               rollingKey = -1;
        for (std::size_t pos = 20; pos <= limit; ++pos) {
            rollingKey = seq2int(seq, static_cast<int>(pos), kKeyLen, rollingKey);
            if (rollingKey >= 0) {
                ++kmerCounts[static_cast<std::size_t>(rollingKey)];
            }
        }
    }

    // Ignore poly‑A (AAAAAAAAAA)
    kmerCounts[0] = 0;

    // Complexity filter
    const auto isLowComplexity = [=](int key) noexcept {
        assert(key >= 0);

        std::array<int, 4> baseCnt = {{0, 0, 0, 0}};
        const auto         ukey    = static_cast<unsigned int>(key);

        for (std::size_t i = 0; i < static_cast<std::size_t>(kKeyLen); ++i) {
            ++baseCnt[(ukey >> (i * 2U)) & 0x3U];
        }

        for (int base : baseCnt) {
            // >= k-4 identical bases
            if (base >= kKeyLen - 4) {
                return true;
            }
        }

        // Too GC-rich or starts with 4xG (GGGG...)
        return (baseCnt[2U] + baseCnt[3U] >= kKeyLen - 2) || ((ukey >> 12U) == 0xFFU);
    };

    // Collect top-N enriched non-trivial k-mers
    struct Candidate {
        int          key   = 0;
        unsigned int count = 0;

        Candidate() = default;
        explicit Candidate(int key_, unsigned int count_)
            : key(key_)
            , count(count_) {}
    };

    std::array<Candidate, kTopCandidates> top {};

    std::uint64_t totalCounts = 0;
    for (std::size_t k = 0; k < kKeySpace; ++k) {
        unsigned int count = kmerCounts[k];
        // We skip trivial k-mers
        if (count == 0 || isLowComplexity(static_cast<int>(k))) {
            continue;
        }
        totalCounts += count;

        // Insert into descending-ordered fixed-size array
        for (int i = kTopCandidates - 1; i >= 0; --i) {
            if (count < top[i].count) {
                // Found where count should go, insert at position+1
                if (i < kTopCandidates - 1) {
                    // Shift elements to make room
                    for (int m = kTopCandidates - 1; m > i + 1; --m) {
                        top[m] = top[m - 1];
                    }
                    top[i + 1] = Candidate {static_cast<int>(k), count};
                }
                break;
            }

            if (i == 0) {
                // count >= all elements, insert at the top
                for (int m = kTopCandidates - 1; m > 0; --m) {
                    top[m] = top[m - 1];
                }
                top[0] = Candidate {static_cast<int>(k), count};
            }
        }
    }

    // Build adapter around each seed
    for (const auto& cand : top) {
        // If the array is prefilled with zeros, we skip entries
        if (cand.key == 0) {
            continue;
        }

        // NOTE: We name this cnt for now to avoid shadowing, this should be moved to another
        // function
        const std::int64_t cnt = cand.count;
        // Check if the candidate is enriched enough, if not we stop looking at smaller ones
        if (cnt < 10
            || cnt * static_cast<std::int64_t>(kKeySpace)
                   < static_cast<std::int64_t>(totalCounts)
                         * static_cast<std::int64_t>(kFoldThreshold)) {
            break;  // not enriched enough
        }

        // NOTE: This copy might be expensive depending on how much this is ran
        const std::string seedSeq = int2seq(static_cast<unsigned int>(cand.key), kKeyLen);

        // Reject seeds with <3 base transitions (low complexity)
        int transitions = std::inner_product(seedSeq.begin(),
                                             seedSeq.end() - 1,
                                             seedSeq.begin() + 1,
                                             0,
                                             std::plus<int>(),
                                             std::not_equal_to<char>());
        if (transitions < 3) {
            continue;
        }

        const auto adapter = getAdapterWithSeed(cand.key, reads, kKeyLen);
        if (!adapter.empty()) {
            return adapter;  // Success, de-novo adapter inferred
        }
    }

    return {};  // No adapter found
}

std::string Evaluator::getAdapterWithSeed(int                                       seed,
                                          const std::vector<std::unique_ptr<Read>>& reads,
                                          int                                       keylen) {
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
    static constexpr std::array<char, 4> kBases{{'A', 'T', 'C', 'G'}};

    // Create a string with length `seqlen`, initialized with 'N' (unknown base)
    std::string result(seqlen, 'N');

    // Convert the integer to a sequence of DNA bases (2 bits per base)
    for (int i = seqlen - 1; i >= 0; --i) {
        // Extract the lowest 2 bits of `val` (mask with 0b11 = 0x3)
        // These 2 bits represent one DNA base (since 2 bits can encode 4 values)
        result[i] = kBases[val & 0x3U]; // TODO: if 0x3U is used elsewhere make it a file level constant
    
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
