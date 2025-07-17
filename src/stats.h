#ifndef STATS_H
#define STATS_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>
#include "read.h"
#include "options.h"


class Stats{
public:
    enum class CurveKey: std::uint8_t { A, T, C, G, N, Mean, GC };

    static constexpr std::size_t InvalidIndex = std::numeric_limits<std::size_t>::max();

    static constexpr std::size_t QualityCurveSize = 5;
    static constexpr std::size_t ContentCurveSize = 6;
    static constexpr std::size_t CurveIndexSize   = 7;

    static auto qualityCurveIndex(CurveKey key) noexcept -> std::size_t;
    static auto contentCurveIndex(CurveKey key) noexcept -> std::size_t;


    // this @guessedCycles parameter should be calculated using the first several records
    Stats(Options* opt, bool isRead2 = false, int guessedCycles = 0, int bufferMargin = 1024);
    ~Stats() = default;

    // Copy and move constructors/assignment operators
    Stats(const Stats& other);
    Stats(Stats&& other) noexcept;
    auto operator=(const Stats& other) -> Stats&;
    auto operator=(Stats&& other) noexcept -> Stats&;

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
    void reportJson(ofstream& ofs, const string& padding);
    // a port of HTML report
    void reportHtml(ofstream& ofs, const string& filteringType, const string& readName);
    void reportHtmlQuality(ofstream& ofs, const string& filteringType, const string& readName);
    void reportHtmlContents(ofstream& ofs, const string& filteringType, const string& readName);
    void reportHtmlKMER(ofstream& ofs, const string& filteringType, const string& readName);
    void reportHtmlORA(ofstream& ofs, const string& filteringType, const string& readName);
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

    std::array<std::vector<double>, QualityCurveSize> mQualityCurves;
    std::array<std::vector<double>, ContentCurveSize> mContentCurves;
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