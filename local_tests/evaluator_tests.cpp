#include <gtest/gtest.h>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <memory>
#include <filesystem>

#define private public
#include "evaluator.h"
#undef private
#include "fastqreader.h"
#include "read.h"
#include "options.h"

using std::string;
using std::unordered_map;
using std::unique_ptr;
using std::vector;

static string create_fastq(const vector<string>& names, const vector<string>& seqs, const string& path) {
    std::ofstream ofs(path);
    for(size_t i=0;i<names.size();++i){
        ofs << names[i] << "\n";
        ofs << seqs[i] << "\n";
        ofs << "+\n";
        ofs << string(seqs[i].size(),'I') << "\n";
    }
    return path;
}

TEST(EvaluatorTests, ComputeSeqLen) {
    vector<string> names = {"@r1","@r2"};
    vector<string> seqs = {"ATCG", "ATCGAT"};
    string path = create_fastq(names, seqs, "len.fq");
    Evaluator eval(nullptr);
    EXPECT_EQ(6, eval.computeSeqLen(path));
}

TEST(EvaluatorTests, EvaluateReadNum) {
    vector<string> names = {"@r1","@r2","@r3"};
    vector<string> seqs(3, "ATCG");
    string path = create_fastq(names, seqs, "readnum.fq");
    Options opt; opt.in1 = path;
    Evaluator eval(&opt);
    long readNum = 0;
    eval.evaluateReadNum(readNum);
    EXPECT_EQ(3, readNum);
}

TEST(EvaluatorTests, MatchKnownAdapter) {
    string adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    EXPECT_EQ(adapter, Evaluator::matchKnownAdapter(adapter+"AAAA"));
    EXPECT_EQ("", Evaluator::matchKnownAdapter("TTTTTTTTTT"));
}

TEST(EvaluatorTests, SeqIntRoundTrip) {
    Evaluator eval(nullptr);
    string s = "ATCGAT";
    int val = eval.seq2int(s,0,static_cast<int>(s.size()),-1);
    EXPECT_EQ(s, eval.int2seq(val, static_cast<int>(s.size())));
}

TEST(EvaluatorTests, Seq2IntInvalidBase) {
    Evaluator eval(nullptr);
    EXPECT_EQ(-1, eval.seq2int("ATCN",0,4,-1));
}