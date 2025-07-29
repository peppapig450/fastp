#include "adaptertrimmer.h"
#include "matcher.h"
#include "read.h"

#include <algorithm>

// Anonymous namespace here allows internal only linkage of these two helpers
namespace {  // (anon)

// REVIEW: these might be better as std::size_t
// One mismatch allowed for every eight bases compared
constexpr int MismatchStride = 8;

inline auto allowedMismatch(int length) noexcept -> int { return length / MismatchStride; }
}  // namespace

// TODO: Taking an overlap analysis object as input is probably better here
auto AdapterTrimmer::trimByOverlapAnalysis(Read*         r1,
                                           Read*         r2,
                                           FilterResult* fr,
                                           int           diffLimit,
                                           int           overlapRequire,
                                           double        diffPercentLimit) -> bool {
    OverlapResult ov =
        OverlapAnalysis::analyze(r1, r2, diffLimit, overlapRequire, diffPercentLimit);

    return trimByOverlapAnalysis(r1, r2, fr, ov);
}

// TODO: look into taking references here instead of pointers
auto AdapterTrimmer::trimByOverlapAnalysis(Read*         r1,
                                           Read*         r2,
                                           FilterResult* fr,
                                           OverlapResult ov,
                                           int           frontTrimmed1,
                                           int           frontTrimmed2) -> bool {
    if (!ov.overlapped || ov.offset >= 0) {
        return false;
    }

    // REVIEW: This is a bit redundant a lambda MIGHT be better here for clarity, but since
    // it is only two variables maybe not.
    const auto r1Len      = static_cast<int>(r1->length());
    const auto r2Len      = static_cast<int>(r2->length());
    const int  overlapLen = ov.overlap_len;

    // 5'      ......frontTrimmed1......|------------------------------------------|----- 3'
    // 3' -----|-------------------------------------------|......frontTrimmed2.....      5'

    const auto trimLenR1 = std::min(r1Len, overlapLen + frontTrimmed2);
    const auto trimLenR2 = std::min(r2Len, overlapLen + frontTrimmed1);

    const std::string adapter1 = r1->seq().substr(trimLenR1);
    const std::string adapter2 = r2->seq().substr(trimLenR2);

#ifdef _DEBUG
    // clang-format off
    std::cerr << adapter1   << '\n'
              << adapter2   << '\n'
              << "frontTrimmed1 = " << frontTrimmed1 << '\n'
              << "frontTrimmed2 = " << frontTrimmed2 << '\n'
              << "overlap = "       << ov.offset     << ", " 
              << overlapLen << ", " << ov.diff       << '\n';
    // clang-format on
    r1->print();
    r2->reverseComplement().print();
    std::cerr << '\n';
#endif  // _DEBUG

    r1->resize(trimLenR1);
    r2->resize(trimLenR2);
    fr->addAdapterTrimmed(adapter1, adapter2);
    return true;
}

auto AdapterTrimmer::trimByMultiSequences(Read*                     r,
                                          FilterResult*             fr,
                                          std::vector<std::string>& adapterList,
                                          bool                      isR2,
                                          bool                      incTrimmedCounter) -> bool {
    if (adapterList.empty()) {
        return false;
    }

    // TODO: this would probably be better as std::size_t
    const auto listSize = static_cast<int>(adapterList.size());
    // clang-format off
    const int matchReq = listSize > 256 ? 6 :
                         listSize >  16 ? 5 : 4;
    // clang-format on

    const auto& originalSeq = r->seq();
    bool        trimmed     = false;

    for (auto& adapter : adapterList) {
        trimmed |= trimBySequence(r, nullptr, adapter, isR2, matchReq);
    }

    if (trimmed) {
        std::string adapter = originalSeq.substr(r->length());
        if (fr != nullptr) {
            fr->addAdapterTrimmed(adapter, isR2, incTrimmedCounter);
        } else {
            std::cerr << adapter << '\n';
        }
    }

    return trimmed;
}

// TODO: matchReq is probably better as std::size_t
auto AdapterTrimmer::trimBySequence(Read*              r,
                                    FilterResult*      fr,
                                    const std::string& adapterSeq,
                                    bool               isR2,
                                    int                matchReq) -> bool {
    const auto adapterLength = static_cast<int>(adapterSeq.size());
    // LATER: figure out what Req means
    if (adapterLength < matchReq) {
        return false;
    }

    // NOTE: if matchReq is std::size_t this should be too
    const auto  readLength  = static_cast<int>(r->length());
    const char* adapterData = adapterSeq.data();
    const char* readData    = r->seq().data();

    // we start from negative numbers since the Illumina adapter dimer usually have the first A
    // skipped as A-tailing
    int startPos = (adapterLength >= 16)   ? -4
                   : (adapterLength >= 12) ? -3
                   : (adapterLength >= 8)  ? -2
                                           : 0;

    bool found = false;
    int  pos   = 0;
    for (pos = startPos; pos < readLength - matchReq && !found; ++pos) {
        const int readShift    = std::max(0, pos);
        const int adapterShift = std::max(0, -pos);
        const int readAvail    = readLength - readShift;
        const int adapterAvail = adapterLength - adapterShift;

        const char* left  = readData + readShift;
        const char* right = adapterData + adapterShift;

        // We let the function see at most one extra base to accommodate a single
        // insertion/deletion
        const int maxCmpLenA = readAvail;
        const int maxCmpLenB = adapterAvail;
        const int cmpLenMin  = std::min(maxCmpLenA, maxCmpLenB);

        const int allowedMismatchCount = allowedMismatch(cmpLenMin);

        const int diff =
            Matcher::editDistance(left, cmpLenMin, right, cmpLenMin, allowedMismatchCount);

        // if we're aligning the *full* adapter length, only allow a perfect match;
        // partial (prefix) matches may still have up to allowedMismatchCount mismatches
        bool fullAdapter = (cmpLenMin == adapterLength);
        if (diff >= 0 && diff <= allowedMismatchCount && (!fullAdapter || diff == 0)) {
            found = true;
            break;
        }
    }

    if (found) {
        if (pos < 0) {
            std::string adapter = adapterSeq.substr(0, adapterLength + pos);
            r->resize(0);
            if (fr != nullptr) {
                fr->addAdapterTrimmed(adapter, isR2);
            }

        } else {
            std::string adapter = r->seq().substr(pos, readLength - pos);
            r->resize(pos);
            if (fr != nullptr) {
                fr->addAdapterTrimmed(adapter, isR2);
            }
        }
        return true;
    }

    return false;
}

bool AdapterTrimmer::test() {
    Read r("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////E");
    std::string adapter = "TTTTCCACGGGGATACTACTG";
    bool trimmed = AdapterTrimmer::trimBySequence(&r, nullptr, adapter);

    const std::string expected1 = "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA";
    if (*r.mSeq != expected1) {
        std::cerr << "[Test 1 Failed]\n"
                  << "Expected: " << expected1 << "\n"
                  << "Actual:   " << *r.mSeq << "\n";
        return false;
    }

    Read read("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGGAAATTTCCCGGGAAATTTCCCGGGATCGATCGATCGATCGAATTCC",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");
    std::vector<std::string> adapterList = {
        "GCTAGCTAGCTAGCTA",
        "AAATTTCCCGGGAAATTTCCCGGG",
        "ATCGATCGATCGATCG",
        "AATTCCGGAATTCCGG"
    };
    trimmed = AdapterTrimmer::trimByMultiSequences(&read, nullptr, adapterList);

    const std::string expected2 = "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG";
    if (*read.mSeq != expected2) {
        std::cerr << "[Test 2 Failed]\n"
                  << "Expected: " << expected2 << "\n"
                  << "Actual:   " << *read.mSeq << "\n";
        return false;
    }

    return true;
}