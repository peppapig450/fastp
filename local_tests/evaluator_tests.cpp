#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#define private public
#include "evaluator.h"
#undef private
#include "adapter_reads.h"
#include "fastqreader.h"
#include "options.h"

using std::string;
using std::unique_ptr;
using std::unordered_map;
using std::vector;

static string create_fastq(const vector<string>& names,
                           const vector<string>& seqs,
                           const string&         path) {
    std::ofstream ofs(path);
    for (size_t i = 0; i < names.size(); ++i) {
        ofs << names[i] << "\n";
        ofs << seqs[i] << "\n";
        ofs << "+\n";
        ofs << string(seqs[i].size(), 'I') << "\n";
    }
    return path;
}

TEST(EvaluatorTests, ComputeSeqLen) {
    vector<string> names = {"@r1", "@r2"};
    vector<string> seqs  = {"ATCG", "ATCGAT"};
    string         path  = create_fastq(names, seqs, "len.fq");
    Evaluator      eval(nullptr);
    EXPECT_EQ(6, eval.computeSeqLen(path));
}

TEST(EvaluatorTests, EvaluateReadNum) {
    vector<string> names = {"@r1", "@r2", "@r3"};
    vector<string> seqs(3, "ATCG");
    string         path = create_fastq(names, seqs, "readnum.fq");
    Options        opt;
    opt.in1 = path;
    Evaluator eval(&opt);
    long      readNum = 0;
    eval.evaluateReadNum(readNum);
    EXPECT_EQ(3, readNum);
}

TEST(EvaluatorTests, MatchKnownAdapter) {
    string adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    EXPECT_EQ(adapter, Evaluator::matchKnownAdapter(adapter + "AAAA"));
    EXPECT_EQ("", Evaluator::matchKnownAdapter("TTTTTTTTTT"));
}

TEST(EvaluatorTests, CheckKnownAdaptersPositive) {
    Evaluator eval(nullptr);
    auto      reads = createReadsWithKnownAdapter();
    EXPECT_EQ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", eval.checkKnownAdapters(reads));
}

TEST(EvaluatorTests, CheckKnownAdaptersNegative) {
    Evaluator eval(nullptr);
    auto      reads = createReadsWithoutAdapter();
    EXPECT_EQ("", eval.checkKnownAdapters(reads));
}

TEST(EvaluatorTests, SeqIntRoundTrip) {
    Evaluator eval(nullptr);
    string    s   = "ATCGAT";
    int       val = eval.seq2int(s, 0, static_cast<int>(s.size()), -1);
    EXPECT_EQ(s, eval.int2seq(val, static_cast<int>(s.size())));
}

TEST(EvaluatorTests, Seq2IntInvalidBase) {
    Evaluator eval(nullptr);
    EXPECT_EQ(-1, eval.seq2int("ATCN", 0, 4, -1));
}

TEST(EvaluatorTests, RollingSeq2IntMatchesFreshCall) {
    Evaluator eval(nullptr);
    string    seq     = "ATCGATCG";
    int       keylen  = 4;
    int       rolling = -1;
    for (int i = 0; i <= static_cast<int>(seq.size()) - keylen; ++i) {
        rolling         = eval.seq2int(seq, i, keylen, rolling);
        int fromScratch = eval.seq2int(seq, i, keylen, -1);
        EXPECT_EQ(fromScratch, rolling) << "Mismatch at position " << i;
    }
}

TEST(EvaluatorTests, IsTwoColorSystem) {
    // FASTQ with read name starting with @NS should be detected as two-color
    vector<string> namesTrue = {"@NS123"};
    vector<string> seqs      = {"ATCG"};
    string         pathTrue  = create_fastq(namesTrue, seqs, "twocolor_true.fq");
    Options        optTrue;
    optTrue.in1 = pathTrue;
    Evaluator evalTrue(&optTrue);
    EXPECT_TRUE(evalTrue.isTwoColorSystem());

    // FASTQ with non-matching read name should not be detected as two-color
    vector<string> namesFalse = {"@XX123"};
    string         pathFalse  = create_fastq(namesFalse, seqs, "twocolor_false.fq");
    Options        optFalse;
    optFalse.in1 = pathFalse;
    Evaluator evalFalse(&optFalse);
    EXPECT_FALSE(evalFalse.isTwoColorSystem());
}

TEST(EvaluatorTests, EvaluateSeqLenAndOverRepSeqs) {
    const int len1   = 12;
    const int len2   = 20;
    const int reads1 = 130;
    const int reads2 = 60;

    vector<string> names1;
    vector<string> seqs1;
    for (int i = 0; i < reads1; ++i) {
        names1.push_back("@A" + std::to_string(i));
        seqs1.emplace_back(len1, 'A');
    }

    vector<string> names2;
    vector<string> seqs2;
    for (int i = 0; i < reads2; ++i) {
        names2.push_back("@C" + std::to_string(i));
        seqs2.emplace_back(len2, 'C');
    }

    string path1 = create_fastq(names1, seqs1, "overrep1.fq");
    string path2 = create_fastq(names2, seqs2, "overrep2.fq");

    Options opt;
    opt.in1 = path1;
    opt.in2 = path2;

    Evaluator eval(&opt);
    eval.evaluateSeqLen();
    EXPECT_EQ(len1, opt.seqLen1);
    EXPECT_EQ(len2, opt.seqLen2);

    eval.evaluateOverRepSeqs();

    auto it1 = opt.overRepSeqs1.find(string(10, 'A'));
    EXPECT_NE(opt.overRepSeqs1.end(), it1);
    EXPECT_GE(it1->second, 500);

    auto it2 = opt.overRepSeqs2.find(string(10, 'C'));
    EXPECT_NE(opt.overRepSeqs2.end(), it2);
    EXPECT_GE(it2->second, 500);
}
