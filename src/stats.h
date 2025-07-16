#ifndef STATS_H
#define STATS_H

#include <array>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>
#include "read.h"
#include "options.h"


class Stats{
public:
    // this @guessedCycles parameter should be calculated using the first several records
    Stats(Options* opt, bool isRead2 = false, int guessedCycles = 0, int bufferMargin = 1024);
    ~Stats() = default;
    int getCycles();
    long getReads();
    long getBases();
    long getQ20();
    long getQ30();
    long getQ40();
    long getGCNumber();
    auto getQualHist() const -> const std::array<long, 128>&;
    // by default the qualified qual score is Q20 ('5')
    void statRead(Read* r);

    static Stats* merge(vector<Stats*>& list);
    void print();
    void summarize(bool forced = false);
    // a port of JSON report
    void reportJson(ofstream& ofs, string padding);
    // a port of HTML report
    void reportHtml(ofstream& ofs, string filteringType, string readName);
    void reportHtmlQuality(ofstream& ofs, string filteringType, string readName);
    void reportHtmlContents(ofstream& ofs, string filteringType, string readName);
    void reportHtmlKMER(ofstream& ofs, string filteringType, string readName);
    void reportHtmlORA(ofstream& ofs, string filteringType, string readName);
    bool isLongRead();
    void initOverRepSeq();
    int getMeanLength();

public:
    static string list2string(double* list, int size);
    static string list2string(double* list, int size, long* coords);
    static string list2string(long* list, int size);
    static int base2val(char base) noexcept;

    // Base types for all count/statistics vectors
    using CountVector = std::vector<long>;

    // Compound types for per-base cycle data
    static constexpr std::size_t BaseCount = 8;
    using BaseArray = std::array<long, BaseCount>;
    using BaseCycleArray = std::array<CountVector, BaseCount>;

private:
    void extendBuffer(int newBufLen);
    string makeKmerTD(int i, int j);
    string kmer3(int val);
    string kmer2(int val);
    void deleteOverRepSeqDist();
    bool overRepPassed(string& seq, long count);

private:
    Options* mOptions;
    bool mIsRead2;
    long mReads;
    int mEvaluatedSeqLen;
    /* 
    why we use 8 here?
    map A/T/C/G/N to 0~7 by their ASCII % 8:
    'A' % 8 = 1
    'T' % 8 = 4
    'C' % 8 = 3
    'G' % 8 = 7
    'N' % 8 = 6
    */
    // TODO: replace this indexing system as it wastes memory
    // TODO: using structs would likely improve cache locality for some of these
    BaseCycleArray mCycleQ30Bases;
    BaseCycleArray mCycleQ20Bases;
    BaseCycleArray mCycleBaseContents;
    BaseCycleArray mCycleBaseQual;

    CountVector mCycleTotalBase;
    CountVector mCycleTotalQual;
    CountVector mKmer;

    // Replace with const?
    std::array<long, 128> mBaseQualHistogram{}; // Initializing at construction like this is preferred

    std::unordered_map<string, std::vector<double>> mQualityCurves;
    std::unordered_map<string, std::vector<double>> mContentCurves;
    std::unordered_map<string, long> mOverRepSeq;
    std::unordered_map<string, std::vector<long>> mOverRepSeqDist;


    int mCycles;
    int mBufLen;
    long mBases;
    BaseArray mQ20Bases;
    BaseArray mQ30Bases;
    BaseArray mBaseContents;
    long mQ20Total;
    long mQ30Total;
    long mQ40Total;
    bool summarized;
    long mKmerMax;
    long mKmerMin;
    int mKmerBufLen;
    long mLengthSum;
};

#endif