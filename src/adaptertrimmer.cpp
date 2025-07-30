#include "adaptertrimmer.h"

#include <algorithm>
#include <cstring>

#include "matcher.h"
#include "read.h"

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

#if _DEBUG
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

auto AdapterTrimmer::trimBySequence(Read*              r,
                                    FilterResult*      fr,
                                    const std::string& adapterSeq,
                                    bool               isR2,
                                    int                matchReq) -> bool {
    const auto adapterLen = static_cast<int>(adapterSeq.size());
    if (adapterLen < matchReq) {
        return false;
    }

    const auto& readSeq     = r->seq();
    const auto  readLen     = static_cast<int>(r->length());
    const char* adapterData = adapterSeq.data();
    const char* readData    = readSeq.data();

    int pos = 0;  // this is filled by locateAdapter
    if (!Matcher::locateAdapter(readData, readLen, adapterData, adapterLen, MismatchStride, matchReq, pos)) {
        return false;
    }

    std::string adapter;
    if (pos < 0) {
        adapter = adapterSeq.substr(0, adapterLen + pos);
        r->resize(0);
    } else {
        adapter = readSeq.substr(pos, readLen - pos);
        r->resize(pos);
    }

    if (fr != nullptr) {
        fr->addAdapterTrimmed(adapter, isR2);
    }
    return true;
}

bool AdapterTrimmer::test() {
    Read r("@name",
           "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG",
           "+",
           "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////E");

    std::string adapter = "TTTTCCACGGGGATACTACTG";
    bool        trimmed = AdapterTrimmer::trimBySequence(&r, nullptr, adapter);

    const std::string expected1 = "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA";
    if (r.seq() != expected1) {
        std::cerr << "[Test 1 Failed]\n"
                  << "Expected: " << expected1 << "\n"
                  << "Actual:   " << *r.mSeq << "\n";
        return false;
    }

    Read read("@name",
              "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGGAAATTTCCCGGGAAATTTCCCGGGATCG"
              "ATCGATCGATCGAATTCC",
              "+",
              "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////"
              "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");

    std::vector<std::string> adapterList = {"GCTAGCTAGCTAGCTA",
                                            "AAATTTCCCGGGAAATTTCCCGGG",
                                            "ATCGATCGATCGATCG",
                                            "AATTCCGGAATTCCGG"};
    trimmed = AdapterTrimmer::trimByMultiSequences(&read, nullptr, adapterList);

    const std::string expected2 = "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG";
    if (read.seq() != expected2) {
        std::cerr << "[Test 2 Failed]\n"
                  << "Expected: " << expected2 << "\n"
                  << "Actual:   " << read.seq() << "\n";
        return false;
    }

    return true;
}