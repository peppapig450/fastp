#pragma once

#include <memory>

#include "read.h"

struct OverlapResult {
    bool overlapped {false};
    // reverse_seq (RC) is shifted right (positive) / left (negative) vs forward_read
    int  offset {0};
    int  overlap_len {0};
    int  diff {0};
    bool hasGap {false};

    // We define our own constructor (ctor) to enable list init of all the contents
    OverlapResult(bool overlapped_, int offset_, int overlap_len_, int diff_, bool has_gap_)
        : overlapped(overlapped_)
        , offset(offset_)
        , overlap_len(overlap_len_)
        , diff(diff_)
        , hasGap(has_gap_) {}

    OverlapResult() = default;
};

class OverlapAnalysis {
public:
    // All members are static so we disallow instantiation
    OverlapAnalysis()  = delete;
    ~OverlapAnalysis() = delete;

    static auto analyze(const std::string& forward_read,
                        const std::string& reverse_rc_read,
                        int                diff_limit,
                        int                min_overlap_len,
                        double             diff_percent_limit,
                        bool               allow_gap = false) -> OverlapResult;

    static auto analyze(const Read& forward_read,
                        const Read& reverse_rc_read,
                        int         diff_limit,
                        int         min_overlap_len,
                        double      diff_percent_limit,
                        bool        allow_gap = false) -> OverlapResult;

    static auto analyze(const Read* forward_read,
                        const Read* reverse_rc_read,
                        int         diff_limit,
                        int         min_overlap_len,
                        double      diff_percent_limit,
                        bool        allow_gap = false) -> OverlapResult;

    static auto merge(const Read& read1, const Read& read2, const OverlapResult& result)
        -> std::unique_ptr<Read>;

    static auto test() -> bool;
};
