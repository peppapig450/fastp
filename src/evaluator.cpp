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

std::string Evaluator::getAdapterWithSeed(int seed, const std::vector<std::unique_ptr<Read>>& loadedReads, int keylen) {
    // we have to shift last cycle for evaluation since it is so noisy, especially for Illumina data
    const int shiftTail = max(1, mOptions->trim.tail1);
    const int MAX_SEARCH_LENGTH = 500;
    NucleotideTree forwardTree(mOptions);
    // forward search
    for(const auto& r : loadedReads) {
        int key = -1;
        const int limit = static_cast<int>(r->length()) - keylen - shiftTail;
        const auto& seq = r->seq();
        for(int pos = 20; pos <= limit && pos <MAX_SEARCH_LENGTH; pos++) {
            key = seq2int(seq, pos, keylen, key);
            if(key == seed) {
                forwardTree.addSeq(seq.substr(pos+keylen, limit - pos));
            }
        }
    }
    bool reachedLeaf = true;
    std::string forwardPath = forwardTree.getDominantPath(reachedLeaf);

    NucleotideTree backwardTree(mOptions);
    // backward search
    for(const auto& r : loadedReads) {
        int key = -1;
        const int limit = static_cast<int>(r->length()) - keylen - shiftTail;
        const auto& seq = r->seq();
        for(int pos = 20; pos <= limit && pos <MAX_SEARCH_LENGTH; pos++) {
            key = seq2int(seq, pos, keylen, key);
            if(key == seed) {
                std::string subseq =  seq.substr(0, pos);
                std::string rcseq = reverse(subseq);
                backwardTree.addSeq(rcseq);
            }
        }
    }
    std::string backwardPath = backwardTree.getDominantPath(reachedLeaf);

    std::string adapter = reverse(backwardPath) + int2seq(seed, keylen) + forwardPath;
    if(adapter.length()>60)
        adapter.resize(60);

    std::string matchedAdapter = matchKnownAdapter(adapter);
    if(!matchedAdapter.empty()) {
        const auto& knownAdapters = adapters::getKnown();
        cerr << knownAdapters.at(matchedAdapter) << endl << matchedAdapter << endl;
        return matchedAdapter;
    } else {
        if(reachedLeaf) {
            cerr << adapter << endl;
            return adapter;
        } else {
            return "";
        }
    }
}

std::string Evaluator::matchKnownAdapter(const std::string& seq) {
    const auto& knownAdapters = adapters::getKnown();
    std::unordered_map<std::string, std::string>::const_iterator iter;
    for(iter = knownAdapters.begin(); iter != knownAdapters.end(); iter++) {
        std::string adapter = iter->first;
        std::string desc = iter->second;
        if(seq.length()<adapter.length()) {
            continue;
        }
        int diff = 0;
        for(int i=0; i<adapter.length() && i<seq.length(); i++) {
            if(adapter[i] != seq[i])
                diff++;
        }
        if(diff == 0)
            return adapter;
    }
    return "";
}

auto Evaluator::int2seq(unsigned int val, int seqlen) const -> std::string {
    char bases[4] = {'A', 'T', 'C', 'G'};
    std::string ret(seqlen, 'N');
    int done = 0;
    while(done < seqlen) {
        ret[seqlen - done - 1] = bases[val & 0x03];
        val = (val >> 2);
        done++;
    }
    return ret;
}

auto Evaluator::seq2int(const std::string& seq, int pos, int keylen, int lastVal) const -> int {
    if (lastVal >= 0) {
        const int mask = (1 << (keylen*2 )) - 1;
        int key = (lastVal<<2) & mask;
        char base = seq[pos + keylen - 1];
        switch (base) {
            case 'A':
                key += 0;
                break;
            case 'T':
                key += 1;
                break;
            case 'C':
                key += 2;
                break;
            case 'G':
                key += 3;
                break;
            default:
                // N or anything else
                return -1;
        }
        return key;
    } else {
        int key = 0;
        for(int i=pos; i<keylen+pos; i++) {
            key = (key << 2);
            char base = seq[i];
            switch (base) {
                case 'A':
                    key += 0;
                    break;
                case 'T':
                    key += 1;
                    break;
                case 'C':
                    key += 2;
                    break;
                case 'G':
                    key += 3;
                    break;
                default:
                    // N or anything else
                    return -1;
            }
        }
        return key;
    }
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
