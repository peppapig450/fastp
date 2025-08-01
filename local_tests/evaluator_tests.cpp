#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <random>
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

static void generate_reads_with_adapter(int                 numReads,
                                        const string&       adapter,
                                        vector<string>&     names,
                                        vector<string>&     seqs) {
    std::mt19937                          rng(1);
    std::uniform_int_distribution<int>    dist(0, 3);
    const char                            bases[4] = {'A', 'T', 'C', 'G'};
    names.clear();
    seqs.clear();
    names.reserve(numReads);
    seqs.reserve(numReads);
    for (int i = 0; i < numReads; ++i) {
        string prefix(20, 'A');
        string suffix(10, 'A');
        for (char& c : prefix) c = bases[dist(rng)];
        for (char& c : suffix) c = bases[dist(rng)];
        names.push_back("@r" + std::to_string(i));
        seqs.push_back(prefix + adapter + suffix);
    }
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

// Ensure evaluateReadNum estimates the total read count using a bytes-per-read
// calculation when the input exceeds the internal sampling limits.
TEST(EvaluatorTests, EvaluateReadNumUsesBytesPerRead) {
    Options opt;
    opt.in1 = MOCK_LARGE_DATASET; // trigger synthetic large dataset
    Evaluator eval(&opt);
    long readNum = 0;
    eval.evaluateReadNum(readNum);

    const long READ_LIMIT = 512 * 1024;
    const long TOTAL_READS = 1024 * 1024; // defined in stub
    const long EXPECTED = static_cast<long>(TOTAL_READS * 1.01);

    EXPECT_GT(readNum, READ_LIMIT);
    EXPECT_NEAR(EXPECTED, readNum, TOTAL_READS * 0.02);
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

TEST(EvaluatorTests, ComputeOverRepSeqPrunesSubstrings) {
    std::mt19937                       rng(1);
    std::uniform_int_distribution<int> dist(0, 3);
    const char                         bases[4] = {'A', 'C', 'G', 'T'};
    auto rand_seq = [&](int len) {
        std::string s;
        s.reserve(len);
        for (int i = 0; i < len; ++i) s += bases[dist(rng)];
        return s;
    };

    std::vector<std::string> names;
    std::vector<std::string> seqs;

    std::string tenC(10, 'C');
    for (int i = 0; i < 500; ++i) {
        names.push_back("@c10_" + std::to_string(i));
        seqs.push_back(tenC + rand_seq(91));
    }

    std::string tenA(10, 'A');
    for (int i = 0; i < 600; ++i) {
        names.push_back("@a10_" + std::to_string(i));
        seqs.push_back(tenA + rand_seq(91));
    }

    std::string twenty = tenA + std::string(10, 'T');
    for (int i = 0; i < 100; ++i) {
        names.push_back("@a20_" + std::to_string(i));
        seqs.push_back(twenty + rand_seq(81));
    }

    std::string forty = rand_seq(40);
    for (int i = 0; i < 20; ++i) {
        names.push_back("@r40_" + std::to_string(i));
        seqs.push_back(forty + rand_seq(61));
    }

    std::string hundred = rand_seq(100);
    for (int i = 0; i < 5; ++i) {
        names.push_back("@r100_" + std::to_string(i));
        seqs.push_back(hundred + rand_seq(1));
    }

    std::string path = create_fastq(names, seqs, "overrep_substring.fq");

    Options    opt;
    opt.in1 = path;
    Evaluator eval(&opt);
    eval.evaluateSeqLen();
    eval.evaluateOverRepSeqs();

    auto& map = opt.overRepSeqs1;

    auto itTenC = map.find(tenC);
    EXPECT_NE(map.end(), itTenC);
    EXPECT_GE(itTenC->second, 500);

    auto itTwenty = map.find(twenty);
    EXPECT_NE(map.end(), itTwenty);
    EXPECT_GE(itTwenty->second, 100);

    auto itForty = map.find(forty);
    EXPECT_NE(map.end(), itForty);
    EXPECT_GE(itForty->second, 20);

    auto itHundred = map.find(hundred);
    EXPECT_NE(map.end(), itHundred);
    EXPECT_GE(itHundred->second, 5);

    EXPECT_EQ(map.end(), map.find(tenA));
}

TEST(EvaluatorTests, EvalAdapterAndReadNumDetectsAdapter) {
    const string adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    vector<string> names;
    vector<string> seqs;
    const int numReads = 10000;
    generate_reads_with_adapter(numReads, adapter, names, seqs);
    string  path = create_fastq(names, seqs, "adapter_eval.fq");
    Options opt;
    opt.in1 = path;
    Evaluator eval(&opt);
    long      readNum = 0;
    string    detected = eval.evalAdapterAndReadNum(readNum, false);
    EXPECT_EQ(adapter, detected);
    EXPECT_NEAR(static_cast<double>(numReads), static_cast<double>(readNum),
                numReads * 0.05);
}

TEST(EvaluatorTests, GetAdapterWithSeedReconstructsAdapter) {
    const string adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    vector<string> names;
    vector<string> seqs;
    const int numReads = 200;
    generate_reads_with_adapter(numReads, adapter, names, seqs);

    vector<unique_ptr<Read>> reads;
    for (int i = 0; i < numReads; ++i) {
        reads.emplace_back(std::make_unique<Read>(names[i], seqs[i]));
    }

    Options    opt;
    Evaluator  eval(&opt);
    const int  keylen = 10;
    int        seed   = eval.seq2int(adapter, 0, keylen, -1);
    string     got    = eval.getAdapterWithSeed(seed, reads, keylen);
    EXPECT_EQ(adapter, got);
}
