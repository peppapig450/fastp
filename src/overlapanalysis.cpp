#include "matcher.h"
#include "overlapanalysis.h"
#include <memory>
#include <string>
#include "sequence.h"

OverlapAnalysis::OverlapAnalysis(){
}


OverlapAnalysis::~OverlapAnalysis(){
}

OverlapResult OverlapAnalysis::analyze(Read* r1, Read* r2, int overlapDiffLimit, int overlapRequire, double diffPercentLimit, bool allowGap) {
    return analyze(r1->mSeq, r2->mSeq, overlapDiffLimit, overlapRequire, diffPercentLimit, allowGap);
}

// ported from the python code of AfterQC
OverlapResult OverlapAnalysis::analyze(string*  r1, string*  r2, int diffLimit, int overlapRequire, double diffPercentLimit, bool allowGap) {
    Sequence rc2 = Sequence(*r2).reverseComplement();
    int len1 = r1->length();
    int len2 = rc2.length();
    // use the pointer directly for speed
    const char* str1 = r1->c_str();
    const char* str2 = rc2.str().c_str();

    int complete_compare_require = 50;

    int overlap_len = 0;
    int offset = 0;
    int diff = 0;

    // forward with no gap
    // a match of less than overlapRequire is considered as unconfident
    while (offset < len1-overlapRequire) {
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1 - offset, len2);
        int overlapDiffLimit = min(diffLimit, (int)(overlap_len * diffPercentLimit));

        diff = 0;
        int i = 0;
        for (i=0; i<overlap_len; i++) {
            if (str1[offset + i] != str2[i]){
                diff += 1;
                if (diff > overlapDiffLimit && i < complete_compare_require)
                    break;
            }
        }
        
        if (diff <= overlapDiffLimit || (diff > overlapDiffLimit && i>complete_compare_require)){
            OverlapResult ov;
            ov.overlapped = true;
            ov.offset = offset;
            ov.overlap_len = overlap_len;
            ov.diff = diff;
            ov.hasGap = false;
            return ov;
        }

        offset += 1;
    }


    // reverse with no gap
    // in this case, the adapter is sequenced since TEMPLATE_LEN < SEQ_LEN
    // check if distance can get smaller if offset goes negative
    // this only happens when insert DNA is shorter than sequencing read length, and some adapter/primer is sequenced but not trimmed cleanly
    // we go reversely
    offset = 0;
    while (offset > -(len2-overlapRequire)){
        // the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1,  len2- abs(offset));
        int overlapDiffLimit = min(diffLimit, (int)(overlap_len * diffPercentLimit));

        diff = 0;
        int i = 0;
        for (i=0; i<overlap_len; i++) {
            if (str1[i] != str2[-offset + i]){
                diff += 1;
                if (diff > overlapDiffLimit && i < complete_compare_require)
                    break;
            }
        }
        
        if (diff <= overlapDiffLimit || (diff > overlapDiffLimit && i>complete_compare_require)){
            OverlapResult ov;
            ov.overlapped = true;
            ov.offset = offset;
            ov.overlap_len = overlap_len;
            ov.diff = diff;
            ov.hasGap = false;
            return ov;
        }

        offset -= 1;
    }

    if(allowGap) {
        // forward with one gap
        offset = 0;
        while (offset < len1-overlapRequire) {
            // the overlap length of r1 & r2 when r2 is move right for offset
            overlap_len = min(len1 - offset, len2);
            int overlapDiffLimit = min(diffLimit, (int)(overlap_len * diffPercentLimit));

            // Ensure the overlap length is valid
            if (overlap_len <= 1 || offset >= len1) {
                offset += 1;
                continue;
            }

            int diff = Matcher::diffWithOneInsertion(str1 + offset, str2, overlap_len-1, overlapDiffLimit);
            if(diff <0 || diff > overlapDiffLimit)
                diff = Matcher::diffWithOneInsertion(str2, str1 + offset, overlap_len-1, overlapDiffLimit);
            
            if (diff <= overlapDiffLimit && diff >=0){
                OverlapResult ov;
                ov.overlapped = true;
                ov.offset = offset;
                ov.overlap_len = overlap_len;
                ov.diff = diff;
                ov.hasGap = true;
                return ov;
            }

            offset += 1;
        }

        // reverse with one gap
        offset = 0;
        while (offset > -(len2-overlapRequire)){
            // the overlap length of r1 & r2 when r2 is move right for offset
            overlap_len = min(len1,  len2- abs(offset));
            int overlapDiffLimit = min(diffLimit, (int)(overlap_len * diffPercentLimit));

            // Ensure we have a valid overlap length and don't exceed bounds
            if (overlap_len <= 1 || std::abs(offset) >= len2) {
                offset -= 1;
                continue;
            }

            // Ensure we do not go beyond string bounds
            const char* shifted_str2 = str2 + std::abs(offset);

            int diff = Matcher::diffWithOneInsertion(str1, shifted_str2, overlap_len-1, overlapDiffLimit);
            if(diff <0 || diff > overlapDiffLimit)
                diff = Matcher::diffWithOneInsertion(shifted_str2, str1, overlap_len-1, overlapDiffLimit);
            
            if (diff <= overlapDiffLimit && diff >=0){
                OverlapResult ov;
                ov.overlapped = true;
                ov.offset = offset;
                ov.overlap_len = overlap_len;
                ov.diff = diff;
                ov.hasGap = true;
                return ov;
            }

            offset -= 1;
        }
    }

    OverlapResult ov;
    ov.overlapped = false;
    ov.offset = ov.overlap_len = ov.diff = 0;
    ov.hasGap = false;
    return ov;
}

Read* OverlapAnalysis::merge(Read* r1, Read* r2, OverlapResult ov) {
    int ol = ov.overlap_len;
    if(!ov.overlapped)
        return NULL;

    int len1 = ol + max(0, ov.offset);
    int len2 = 0; 
    if(ov.offset > 0)
        len2 = r2->length() - ol;

    Read* rr2 = r2->reverseComplement();
    string mergedSeq = r1->mSeq->substr(0, len1);
    if(ov.offset > 0) {
        mergedSeq += rr2->mSeq->substr(ol, len2);
    }

    string mergedQual = r1->mQuality->substr(0, len1);
    if(ov.offset > 0) {
        mergedQual += rr2->mQuality->substr(ol, len2);
    }

    delete rr2;

    string name = *(r1->mName) + " merged_" + to_string(len1) + "_" + to_string(len2);
    string strand = *(r1->mStrand);
    if (strand != "+") {
      strand = strand + " merged_" + to_string(len1) + "_" + to_string(len2);
    }
    Read* mergedRead = new Read(new string(name), new string(mergedSeq), new string(strand), new string(mergedQual));

    return mergedRead;
}

bool OverlapAnalysis::test() {
    // Sequence
    // r1("CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGCCGCTGGAGGTCTCCC");
    // Sequence
    // r2("ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGCCCGTAGGCGCGGCTCCC");

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
        auto s1 = std::unique_ptr<std::string>(new std::string(tcase.seq1));
        auto s2 = std::unique_ptr<std::string>(new std::string(tcase.seq2));
        auto q1 = std::unique_ptr<std::string>(new std::string(tcase.qual1));
        auto q2 = std::unique_ptr<std::string>(new std::string(tcase.qual2));

        // run analyze
        auto overlapResult = OverlapAnalysis::analyze(s1.get(),
                                                      s2.get(),
                                                      tcase.min_ol,
                                                      tcase.match_thr,
                                                      tcase.diff_thr,
                                                      tcase.allow_gap);

        bool resultMatchExpectations =
            overlapResult.overlapped == tcase.overlapped && overlapResult.offset == tcase.offset
            && overlapResult.overlap_len == tcase.overlap_len && overlapResult.diff == tcase.diff
            && overlapResult.hasGap == tcase.has_gap;

        // Wrap into Reads and merge
        Read r1 {new std::string("n1"), s1.release(), new std::string("+"), q1.release()};
        Read r2 {new std::string("n2"), s2.release(), new std::string("+"), q2.release()};
        std::unique_ptr<Read> merged {OverlapAnalysis::merge(&r1, &r2, overlapResult)};

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
