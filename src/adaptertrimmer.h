#pragma once

#include <string>

#include "filterresult.h"
#include "overlapanalysis.h"

class Read;  // Forward declaration

// NOTE: This can probably be a namespace, or can use members for reduced boilerplate
class AdapterTrimmer {
public:
    AdapterTrimmer()  = default;
    ~AdapterTrimmer() = default;

    static auto trimByOverlapAnalysis(Read*         r1,
                                      Read*         r2,
                                      FilterResult* fr,
                                      int           diffLimit,
                                      int           overlapRequire,
                                      double        diffPercentLimit) -> bool;
    static auto trimByOverlapAnalysis(Read*         r1,
                                      Read*         r2,
                                      FilterResult* fr,
                                      OverlapResult ov,
                                      int           frontTrimmed1 = 0,
                                      int           frontTrimmed2 = 0) -> bool;
    static auto trimBySequence(Read*              r,
                               FilterResult*      fr,
                               const std::string& adapterSeq,
                               bool               isR2     = false,
                               int                matchReq = 4) -> bool;
    static auto trimByMultiSequences(Read*                     r1,
                                     FilterResult*             fr,
                                     std::vector<std::string>& adapterList,
                                     bool                      isR2              = false,
                                     bool                      incTrimmedCounter = true) -> bool;
    static auto test() -> bool;
};