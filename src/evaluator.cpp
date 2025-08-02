#include "evaluator.h"

#include <algorithm>
#include <cstdint>
#include <cstddef>
#include <map>
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
    FastqReader reader(mOptions->in1);

    const long READ_LIMIT = 512*1024;
    const long BASE_LIMIT = 151 * 512*1024;
    long records = 0;
    long bases = 0;
    std::size_t firstReadPos = 0;

    std::size_t bytesRead;
    std::size_t bytesTotal;

    bool reachedEOF = false;
    bool first = true;
    while(records < READ_LIMIT && bases < BASE_LIMIT) {
        Read* r = reader.read();
        if(!r) {
            reachedEOF = true;
            break;
        }
        if(first) {
            reader.getBytes(bytesRead, bytesTotal);
            firstReadPos = bytesRead;
            first = false;
        }
        records++;
        bases += r->length();
        delete r;
    }

    readNum = 0;
    if(reachedEOF){
        readNum = records;
    } else if(records>0) {
        // by the way, update readNum so we don't need to evaluate it if splitting output is enabled
        reader.getBytes(bytesRead, bytesTotal);
        double bytesPerRead = (double)(bytesRead - firstReadPos) / (double) records;
        // increase it by 1% since the evaluation is usually a bit lower due to bad quality causes lower compression rate
        readNum = (long) (bytesTotal*1.01 / bytesPerRead);
    }
}

std::string Evaluator::checkKnownAdapters(const std::vector<std::unique_ptr<Read>>& reads) {
    std::map<std::string, std::string> knownAdapters = getKnownAdapter();
    std::map<std::string, int> possibleCounts;
    std::map<std::string, int> mismatches;

    // TODO: These should probably be constexpr
    // for performance, up to 100k reads and 100M bases
    const int MAX_CHECK_READS = 100000;
    const int MAX_CHECK_BASES = MAX_CHECK_READS*1000;
    // if we hit for 1000 times, then we exit
    const int MAX_HIT = 1000;

    const int matchReq = 8;
    const int allowOneMismatchForEach = 16;

    std::map<std::string, std::string>::iterator iter;

    for(iter = knownAdapters.begin(); iter!= knownAdapters.end(); iter++) {
        std::string adapter = iter->first;
        possibleCounts[adapter] = 0;
        mismatches[adapter] = 0;
    }

    long checkedReads = 0;
    long checkedBases = 0;
    int curMaxCount = 0;
    for (const auto& r : reads) {
        const char* rdata = r->seq().data();
        int rlen = r->length();

        ++checkedReads;
        checkedBases += rlen;
        if (checkedReads > MAX_CHECK_READS || checkedBases > MAX_CHECK_BASES || curMaxCount > MAX_HIT) {
            break;
        }

        for(iter = knownAdapters.begin(); iter!= knownAdapters.end(); iter++) {
            std::string adapter = iter->first;
            const char* adata = adapter.c_str();
            int alen = adapter.length();

            if(alen >= rlen) {
                continue;
            }
            // this one is not the candidate, skip it for speedup
            if(curMaxCount > 20 && possibleCounts[adapter] <curMaxCount/10) {

                continue;
            }

            for(int pos = 0; pos<rlen-matchReq; pos++) {
                int cmplen = std::min(rlen - pos, alen);
                int allowedMismatch = cmplen/allowOneMismatchForEach;
                int mismatch = 0;
                bool matched = true;
                for(int i=0; i<cmplen; i++) {
                    if( adata[i] != rdata[i+pos] ){
                        mismatch++;
                        if(mismatch > allowedMismatch) {
                            matched = false;
                            break;
                        }
                    }
                }
                if(matched) {
                    possibleCounts[adapter]++;
                    curMaxCount = std::max(curMaxCount, possibleCounts[adapter]);
                    mismatches[adapter] += mismatch;
                    break;
                }
            }
        }
    }

    std::string adapter;
    int maxCount = 0;
    std::map<std::string, int>::iterator iter2;
    for(iter2 = possibleCounts.begin(); iter2 != possibleCounts.end(); iter2++) {
        if(iter2->second > maxCount) {
            adapter = iter2->first;
            maxCount = iter2->second;
        }
    }
    if(maxCount > checkedReads/50 || (maxCount > checkedReads/200 && mismatches[adapter] < checkedReads)) {
        cerr << knownAdapters[adapter] << endl;
        cerr << adapter << endl;
        return adapter;
    }
    return "";
}

std::string Evaluator::evalAdapterAndReadNum(long& readNum, bool isR2) {
    std::string filename = mOptions->in1;
    if(isR2)
        filename = mOptions->in2;
    FastqReader reader(filename);
    // stat up to 256K reads
    const long READ_LIMIT = 256*1024;
    const long BASE_LIMIT = 151 * READ_LIMIT;
    long records = 0;
    long bases = 0;
    std::size_t firstReadPos = 0;

    std::size_t bytesRead;
    std::size_t bytesTotal;

    std::vector<std::unique_ptr<Read>> loadedReads;
    loadedReads.reserve(READ_LIMIT);
    bool reachedEOF = false;
    bool first = true;

    while(records < READ_LIMIT && bases < BASE_LIMIT) {
        std::unique_ptr<Read> r{reader.read()}; // TODO: make reader return unique_ptr directly
        if(!r) {
            reachedEOF = true;
            break;
        }
        if(first) {
            reader.getBytes(bytesRead, bytesTotal);
            firstReadPos = bytesRead;
            first = false;
        }
        int rlen = r->length();
        bases += rlen;
        loadedReads.emplace_back(std::move(r));
        records++;
    }
    loadedReads.shrink_to_fit();

    readNum = 0;
    if(reachedEOF){
        readNum = records;
    } else if(records>0) {
        // by the way, update readNum so we don't need to evaluate it if splitting output is enabled
        reader.getBytes(bytesRead, bytesTotal);
        double bytesPerRead = (double)(bytesRead - firstReadPos) / (double) records;
        // increase it by 1% since the evaluation is usually a bit lower due to bad quality causes lower compression rate
        readNum = (long) (bytesTotal*1.01 / bytesPerRead);
    }

    // we need at least 10000 valid records to evaluate
    if(records < 10000) {
        return "";
    }

    // TODO: probably want to pass by reference to mutate in place and avoid copying
    std::string knownAdapter = checkKnownAdapters(loadedReads);
     if(knownAdapter.size() > 8) {
        return knownAdapter;
    }

    // we have to shift last cycle for evaluation since it is so noisy, especially for Illumina data
    const std::size_t shiftTail = std::max<std::size_t>(1, mOptions->trim.tail1);

    // why we add trim_tail here? since the last cycle are usually with low quality and should be trimmed
    const int keylen = 10;
    const std::size_t size = 1U << (static_cast<std::size_t>(keylen) * 2);
    std::vector<unsigned int> counts(size, 0);

    for (const auto& r : loadedReads) {
        int key = -1;
        const std::size_t limit = r->length() - static_cast<std::size_t>(keylen) - shiftTail;
        for(std::size_t pos = 20; pos <= limit; pos++) {
            key = seq2int(r->seq(), static_cast<int>(pos), static_cast<int>(keylen), key);
            if(key >= 0) {
                counts[key]++;
            }
        }
    }

    // set AAAAAAAAAA = 0;
    counts[0] = 0;

    // get the top N
    const int topnum = 10;
    int topkeys[topnum] = {0};
    long total = 0;
    for(int k=0; k<size; k++) {
        int atcg[4] = {0};
        for(int i=0; i<keylen; i++) {
            int baseOfBit = (k >> (i*2)) & 0x03;
            atcg[baseOfBit]++;
        }
        bool lowComplexity = false;
        for(int b=0; b<4; b++) {
            if(atcg[b] >= keylen-4)
                lowComplexity=true;
        }
        if(lowComplexity)
            continue;
        // too many GC
        if(atcg[2] + atcg[3] >= keylen-2)
            continue;

        // starts with GGGG
        if( k>>12 == 0xff)
            continue;

        unsigned int val = counts[k];
        total += val;
        for(int t=topnum-1; t>=0; t--) {
            // reach the middle
            if(val < counts[topkeys[t]]){
                if(t<topnum-1) {
                    for(int m=topnum-1; m>t+1; m--) {
                        topkeys[m] = topkeys[m-1];
                    }
                    topkeys[t+1] = k;
                }
                break;
            } else if(t == 0) { // reach the top
                for(int m=topnum-1; m>t; m--) {
                    topkeys[m] = topkeys[m-1];
                }
                topkeys[t] = k;
            }
        }
    }

    const int FOLD_THRESHOLD = 20;
    for(int t=0; t<topnum; t++) {
        int key = topkeys[t];
        std::string seq = int2seq(key, keylen);
        if(key == 0)
            continue;
        long count = counts[key];
        if(count<10 || count*size < total * FOLD_THRESHOLD)
            break;
        // skip low complexity seq
        int diff = 0;
        for(int s=0; s<seq.length() - 1; s++) {
            if(seq[s] != seq[s+1])
                diff++;
        }
        if(diff <3){
            continue;
        }
        std::string adapter = getAdapterWithSeed(key, loadedReads, keylen);
        if(!adapter.empty()){
            return adapter;
        }
    }

    return "";
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
        map<std::string, std::string> knownAdapters = getKnownAdapter();
        cerr << knownAdapters[matchedAdapter] << endl << matchedAdapter << endl;
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
    map<std::string, std::string> knownAdapters = getKnownAdapter();
    map<std::string, std::string>::iterator iter;
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
