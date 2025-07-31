#pragma once

#include <memory>
#include <string>
#include "options.h"
#include "read.h" 
#include <unordered_map>
#include <vector>

class Evaluator{
public:
    explicit Evaluator(Options* opt);
    ~Evaluator();
    // evaluate how many reads are stored in the input file
    void evaluateReadNum(long& readNum);
    std::string evalAdapterAndReadNum(long& readNum, bool isR2);
    bool isTwoColorSystem();
    void evaluateSeqLen();
    void evaluateOverRepSeqs();
    void computeOverRepSeq(std::string filename, std::unordered_map<std::string, long>& hotseqs, int seqLen);
    int computeSeqLen(std::string filename);

    static bool test();
    static std::string matchKnownAdapter(std::string seq);
private:
    Options* mOptions;
    std::string int2seq(unsigned int val, int seqlen);
    int seq2int(const std::string& seq, int pos, int keylen, int lastVal = -1);
    std::string getAdapterWithSeed(int seed, const std::vector<std::unique_ptr<Read>>& loadedReads, int keylen);
    std::string checkKnownAdapters(const std::vector<std::unique_ptr<Read>>& reads);
};
