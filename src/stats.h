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

    enum class BaseIndex : std::uint8_t { A = 0, T, C, G, N, Count };
    static constexpr std::size_t BaseCount = static_cast<std::size_t>(BaseIndex::Count);

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

    // Static lookup array for base values
    static constexpr std::array<signed char, 256> BaseValueLookup = {
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0,  -1, 2,  -1, -1, -1, 3,  -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, 1,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    };
    static auto base2val(char base) noexcept -> std::int8_t;

    // Base types for all count/statistics vectors
    using CountVector = std::vector<long>;

    // helper container for per-base per-cycle statistics
    struct BaseCycleArray {
      private:
        std::vector<long> data;
        std::size_t bases{0};
        std::size_t bufLen{0};

        auto index(std::size_t cycle, std::size_t base) const noexcept -> std::size_t {
            return (cycle * bases) + base;
        }

      public:
        BaseCycleArray() = default;
        BaseCycleArray(std::size_t base, std::size_t len) : data(len * base, 0), bases(base), bufLen(len) {}

        auto operator()(std::size_t cycle, std::size_t base) -> long& { return data[index(cycle, base)]; }

        auto operator()(std::size_t cycle, std::size_t base) const -> const long& {
            return data[index(cycle, base)];
        }

        void resize(std::size_t base, std::size_t newLen) {
            std::vector<long> newData(newLen * base, 0);
            std::size_t minBase = std::min(base, bases);
            std::size_t minLen  = std::min(newLen, bufLen);

            for (std::size_t cycleIdx = 0; cycleIdx < minLen; ++cycleIdx) {
                for (std::size_t baseIdx = 0; baseIdx < minBase; ++baseIdx) {
                    newData[(cycleIdx * base) + baseIdx] = (*this)(cycleIdx, baseIdx);
                }
            }

            data.swap(newData);
            bases = base;
            bufLen = newLen;
        }

        auto cycleLength() const -> std::size_t { return bufLen; }
        auto baseCount() const -> std::size_t { return bases; }

        auto values() -> std::vector<long>& { return data; }
        auto values() const -> const std::vector<long>& { return data; }
    };

    // Compound types for per-base cycle data
    static constexpr std::array<unsigned char, 256> BaseIndexLookup = {
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 5, 2, 5, 5, 5, 3, 5, 5,
        5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    };
    static auto baseIndex(char base) noexcept -> std::size_t;
    using BaseArray = std::array<long, BaseCount>;

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

    // Per-base statistics
    BaseCycleArray mCycleQ30Bases;
    BaseCycleArray mCycleQ20Bases;
    BaseCycleArray mCycleBaseContents;
    BaseCycleArray mCycleBaseQual;

    CountVector mCycleTotalBase;
    CountVector mCycleTotalQual;
    CountVector mKmer;

    //TODO: Replace with const?
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
