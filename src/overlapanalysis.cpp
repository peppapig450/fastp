#include "overlapanalysis.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>
#include <utility>

#include "cassert"
#include "matcher.h"
#include "read.h"
#include "sequence.h"

// Use anonymous namespace for internal linkage, this avoids static and is more preferred
// in modern C++.
namespace {  // (anon)

// Static scan configuration
struct ScanConfig {
    const char* forward_read;
    const char* reverse_rc_read;
    std::size_t forward_read_len;
    std::size_t reverse_rc_read_len;
    int         min_overlap_len;
    int         diff_limit;
    double      diff_percent_limit;
    bool        reverse; /* false = forward, true = reverse */
};

// Dynamic data for *one* offset
struct IterContext {
    int         offset;
    std::size_t overlap_len;
    int         overlap_diff_limit;
    int         direction_step;
};

// Generic overlap-scanner
// Scans for overlaps using a strategy-based template approach.
// The scanning logic is decoupled from the overlap evaluation using lambdas.
//
// Benefits:
// - Separation of concerns: iteration and evaluation are independent
// - Reusability: shared logic for no-gap and gap-aware variants
// - Zero-cost abstraction via inlined lambdas
// - Compile-time safety through template-based evaluator signatures
// - Easier maintenance: scanning logic centralized
//
// The evaluator lambda takes an IterContext (current offset info) and returns
// true if an overlap is found (to stop early), false to continue scanning.
template <class OffsetEvaluator>
inline auto scanOverlaps(const ScanConfig& config,
                         OverlapResult&    result,
                         OffsetEvaluator&& evaluator) -> bool {
    const int direction_step = config.reverse ? -1 : 1;
    int       offset         = 0;

    const auto min_overlap_len = config.min_overlap_len;
    while (
        (config.reverse && offset > -static_cast<int>(config.reverse_rc_read_len - min_overlap_len))
        || (!config.reverse
            && offset < static_cast<int>(config.forward_read_len - min_overlap_len))) {
        std::size_t overlap_len =
            config.reverse ? std::min<std::size_t>(config.forward_read_len,
                                                   config.reverse_rc_read_len
                                                       - static_cast<std::size_t>(std::abs(offset)))
                           : std::min<std::size_t>(config.forward_read_len - offset,
                                                   config.reverse_rc_read_len);

        const auto overlap_fraction = static_cast<double>(overlap_len) * config.diff_percent_limit;
        const auto overlap_diff_limit =
            std::min(config.diff_limit, static_cast<int>(overlap_fraction));

        // clang-format off
        IterContext context {
            offset,
             overlap_len,
             overlap_diff_limit,
             direction_step
        };
        // clang-format on

        if (evaluator(context, result)) {
            return true;
        }

        offset += direction_step;
    }
    return false;
}

// When offset >= 0: read2_RC is shifted right; offset < 0: shifted left
// Returns true and fill |res| on success, otherwise false.
// TODO: maybe use const string references here
inline auto forwardNoGap(const char*    forward_read,
                         const char*    reverse_rc_read,
                         std::size_t    forward_read_len,
                         std::size_t    reverse_rc_read_len,
                         int            min_overlap_len,
                         int            diff_limit,
                         double         diff_percent_limit,
                         OverlapResult& result,
                         bool           scan_reverse /* false = forward, true = reverse*/) -> bool {
    constexpr int min_full_comparison_len = 50;

    ScanConfig config {forward_read,
                       reverse_rc_read,
                       forward_read_len,
                       reverse_rc_read_len,
                       min_overlap_len,
                       diff_limit,
                       diff_percent_limit,
                       scan_reverse};

    return scanOverlaps(config, result, [=](const IterContext& ctx, OverlapResult& result) -> bool {
        const auto& offset             = ctx.offset;
        const auto& overlap_diff_limit = ctx.overlap_diff_limit;

        int         diff = 0;
        std::size_t i    = 0;  // This allows us to use i outside the loop.
        for (; i < ctx.overlap_len; ++i) {
            const char forward_base = scan_reverse ? forward_read[i] : forward_read[offset + i];
            const char reverse_base =
                scan_reverse ? reverse_rc_read[std::abs(offset) + i] : reverse_rc_read[i];

            if (forward_base != reverse_base) {
                ++diff;
                if (diff > overlap_diff_limit && static_cast<int>(i) < min_full_comparison_len) {
                    break;
                }
            }
        }

        bool accept =
            diff <= overlap_diff_limit
            || (diff > overlap_diff_limit && static_cast<int>(i) > min_full_comparison_len);

        if (accept) {
            // clang-format off
            result = {
                true,
                offset,
                static_cast<int>(ctx.overlap_len),
                diff, false
            };
            // clang-format on
        }
        return accept;
    });
}

// Gap-aware branch using Matcher helpers
inline auto gapPass(const char*    forward_read,
                    const char*    reverse_rc_read,
                    std::size_t    forward_read_len,
                    std::size_t    reverse_rc_read_len,
                    int            min_overlap_len,
                    int            diff_limit,
                    double         diff_percent_limit,
                    OverlapResult& result,
                    bool           scan_reverse) -> bool {
    ScanConfig config {forward_read,
                       reverse_rc_read,
                       forward_read_len,
                       reverse_rc_read_len,
                       min_overlap_len,
                       diff_limit,
                       diff_percent_limit,
                       scan_reverse};

    return scanOverlaps(config, result, [=](const IterContext ctx, OverlapResult& result) -> bool {
        // Early exit if length is less than or equal to one
        if (ctx.overlap_len <= 1) {
            return false;
        }

        const auto& offset             = ctx.offset;
        const auto& overlap_len        = ctx.overlap_len;
        const auto& overlap_diff_limit = ctx.overlap_diff_limit;

        // TODO: investigate whether using constant string references here offers performance gains,
        // that would require updating Matcher to take references instead of pointers, and also
        // eliminates the need to use __restrict for pointer aliasing there.
        const char* left_read = scan_reverse ? forward_read : forward_read + offset;
        const char* right_read =
            scan_reverse ? reverse_rc_read + std::abs(offset) : reverse_rc_read;

        const auto overlap_compare_len = static_cast<int>(overlap_len - 1);

        // TODO: probably should adjust the function params in Matcher::diffWithOneInsertion
        // our usage of left/right here does't really match the 'insertion/normal Data' that it
        // uses.
        int diff = Matcher::diffWithOneInsertion(left_read,
                                                 right_read,
                                                 overlap_compare_len,
                                                 overlap_diff_limit);

        if (diff < 0 || diff > overlap_diff_limit) {
            diff = Matcher::diffWithOneInsertion(right_read,
                                                 left_read,
                                                 overlap_compare_len,
                                                 overlap_diff_limit);
        }

        if (diff >= 0 && diff <= overlap_diff_limit) {
            // clang-format off
            result = {
                true,
                offset,
                static_cast<int>(overlap_len),
                diff,
                true
            };
            // clang-format on
        }
        return false;
    });
}
}  // namespace

auto OverlapAnalysis::analyze(const std::string& forward_read,
                              const std::string& reverse_rc_read,
                              int                diff_limit,
                              int                min_overlap_len,
                              double             diff_percent_limit,
                              bool               allow_gap) -> OverlapResult {
    assert(min_overlap_len > 0 && "min_overlap_len must be positive");

    Sequence    reverse_rc = Sequence(reverse_rc_read).reverseComplement();
    const auto& rev_str    = reverse_rc.str();

    // TODO: figure out if raw pointers here are actually faster
    const char*       forward_str        = forward_read.c_str();
    const char*       reverse_rc_str     = rev_str.c_str();
    const std::size_t forward_str_len    = forward_read.length();
    const std::size_t reverse_rc_str_len = reverse_rc.length();

    OverlapResult result;

    // Forward, no gap
    if (forwardNoGap(forward_str,
                     reverse_rc_str,
                     forward_str_len,
                     reverse_rc_str_len,
                     min_overlap_len,
                     diff_limit,
                     diff_percent_limit,
                     result,
                     /*scan_reverse=*/false)) {
        return result;
    }

    // Reverse, no gap
    if (forwardNoGap(forward_str,
                     reverse_rc_str,
                     forward_str_len,
                     reverse_rc_str_len,
                     min_overlap_len,
                     diff_limit,
                     diff_percent_limit,
                     result,
                     /*scan_reverse=*/true)) {
        return result;
    }

    if (allow_gap) {
        // Forward, single gap
        if (gapPass(forward_str,
                    reverse_rc_str,
                    forward_str_len,
                    reverse_rc_str_len,
                    min_overlap_len,
                    diff_limit,
                    diff_percent_limit,
                    result,
                    /*scan_reverse=*/false)) {
            return result;
        }

        // Reverse, single gap
        if (gapPass(forward_str,
                    reverse_rc_str,
                    forward_str_len,
                    reverse_rc_str_len,
                    min_overlap_len,
                    diff_limit,
                    diff_percent_limit,
                    result,
                    /*scan_reverse=*/true)) {
            return result;
        }
    }

    // no overlap found, we return the default constructed OverlapResult
    return result;
}

// Thin wrapper that adapts Read -> string
auto OverlapAnalysis::analyze(const Read& forward_read,
                              const Read& reverse_rc_read,
                              int         diff_limit,
                              int         min_overlap_len,
                              double      diff_percent_limit,
                              bool        allow_gap) -> OverlapResult {
    return analyze(forward_read.seq(),
                   reverse_rc_read.seq(),
                   diff_limit,
                   min_overlap_len,
                   diff_percent_limit,
                   allow_gap);
}

// Thin wrapper for legacy pointer usage
auto OverlapAnalysis::analyze(const Read* forward_read,
                              const Read* reverse_rc_read,
                              int         diff_limit,
                              int         min_overlap_len,
                              double      diff_percent_limit,
                              bool        allow_gap) -> OverlapResult {
    assert(forward_read && reverse_rc_read);

    return analyze(*forward_read,
                   *reverse_rc_read,
                   diff_limit,
                   min_overlap_len,
                   diff_percent_limit,
                   allow_gap);
}

auto OverlapAnalysis::merge(const Read& read1, const Read& read2, const OverlapResult& result)
    -> std::unique_ptr<Read> {
    if (!result.overlapped) {
        return nullptr;
    }

    const auto overlap_len    = result.overlap_len;
    const auto overlap_offset = result.offset;
    const auto len1           = overlap_len + std::max(0, overlap_offset);
    const auto len2 = (overlap_offset > 0) ? static_cast<int>(read2.length()) - overlap_len : 0;

    // reverse complement read2 for merged tail, then discard
    auto read2_reverse_compl =
        std::unique_ptr<Read>(new Read(std::move(read2.reverseComplement())));

    const auto len1_size = static_cast<std::size_t>(len1);

    auto mergedSeq  = read1.seq().substr(0, len1_size);
    auto mergedQual = read1.quality().substr(0, len1_size);

    const auto overlap_len_size = static_cast<std::size_t>(overlap_len);
    const auto len2_size        = static_cast<std::size_t>(len2);
    if (overlap_offset > 0) {
        mergedSeq += read2_reverse_compl->seq().substr(overlap_len_size, len2_size);
        mergedQual += read2_reverse_compl->quality().substr(overlap_len_size, len2_size);
    }

    auto name   = read1.name() + " merged_" + std::to_string(len1) + "_" + std::to_string(len2);
    auto strand = read1.strand();
    if (strand != "+") {
        strand += "merged_" + std::to_string(len1) + "_" + std::to_string(len2);
    }

    return std::unique_ptr<Read>(new Read(name, mergedSeq, strand, mergedQual));
}

auto OverlapAnalysis::test() -> bool {
    struct TestCase {
        std::string seq1;
        std::string seq2;
        std::string qual1;
        std::string qual2;

        int    min_ol;
        int    match_thr;
        double diff_thr;
        bool   allow_gap;

        // What we expect
        bool        overlapped;
        int         offset;
        int         overlap_len;
        int         diff;
        bool        has_gap;
        std::string merged_seq;
        std::string merged_qual;
    };

    auto runTestCase = [&](TestCase const& tcase) -> bool {
        // build inputs
        const std::string& s1 = tcase.seq1;
        const std::string& s2 = tcase.seq2;
        const std::string& q1 = tcase.qual1;
        const std::string& q2 = tcase.qual2;

        // run analyze
        auto overlapResult = OverlapAnalysis::analyze(s1,
                                                      s2,
                                                      tcase.min_ol,
                                                      tcase.match_thr,
                                                      tcase.diff_thr,
                                                      tcase.allow_gap);

        bool resultMatchExpectations =
            overlapResult.overlapped == tcase.overlapped && overlapResult.offset == tcase.offset
            && overlapResult.overlap_len == tcase.overlap_len && overlapResult.diff == tcase.diff
            && overlapResult.hasGap == tcase.has_gap;

        // Wrap into Reads and merge
        Read                  r1 {"n1", s1, "+", q1};
        Read                  r2 {"n2", s2, "+", q2};
        std::unique_ptr<Read> merged {OverlapAnalysis::merge(r1, r2, overlapResult)};

        if (tcase.overlapped) {
            resultMatchExpectations &= merged && *merged->mSeq == tcase.merged_seq
                                       && *merged->mQuality == tcase.merged_qual;
        } else {
            resultMatchExpectations &= (merged == nullptr);
        }
        return resultMatchExpectations;
    };

    std::array<TestCase, 5> cases {{
        // perfect rc overlap
        // clang-format off
        { "ACGTACGT","ACGTACGT", "HHHHHHHH","IIIIIIII",
          0,4,0.1,false,
          true,  0,8,0,false,
          "ACGTACGT","HHHHHHHH"
        },
        // positive offset
        { "AAAAGGGG","AACCCCTT", "HHHHHHHH","IIIIIIII",
          0,4,0.1,false,
          true,  2,6,0,false,
          "AAAAGGGGTT","HHHHHHHHII"
        },
        // negative offset
        { "GGGGAAAA","TTTTCCCCAA","JJJJJJJJ","KKKKKKKKKK",
          0,4,0.1,false,
          true, -2,8,0,false,
          "GGGGAAAA","JJJJJJJJ"
        },
        // gap allowed
        { "ACGTAACGT","ACGTACGT","LLLLLLLLL","MMMMMMMM",
          1,4,0.3,true,
          true,  0,8,0,true,
          "ACGTAACG","LLLLLLLL"
        },
        // no overlap
        { "AAAA","CCCC","NNNN","OOOO",
          0,3,0.1,false,
          false,0,0,0,false,
          "",""
        }  // clang-format on
    }};

    bool passedTests = true;
    for (auto const& tcase : cases) {
        passedTests &= runTestCase(tcase);
    }

    return passedTests;
}
