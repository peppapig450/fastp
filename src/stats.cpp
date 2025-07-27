#include "stats.h"
#include <memory.h>
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include "util.h"


// TODO: a lot of these functions are not used elsewhere so they could be removed

static const std::array<Stats::CurveKey, Stats::QualityCurveSize> QualityOrder = {Stats::CurveKey::A,
                                                                                  Stats::CurveKey::T,
                                                                                  Stats::CurveKey::C,
                                                                                  Stats::CurveKey::G,
                                                                                  Stats::CurveKey::Mean};

static const std::array<Stats::CurveKey, Stats::ContentCurveSize> ContentOrder = {Stats::CurveKey::A,
                                                                                  Stats::CurveKey::T,
                                                                                  Stats::CurveKey::C,
                                                                                  Stats::CurveKey::G,
                                                                                  Stats::CurveKey::GC};

static const std::array<const char*, Stats::QualityCurveSize> QualityNames = {"A", "T", "C", "G", "mean"};

static const std::array<const char*, Stats::ContentCurveSize> ContentNames = {"A", "T", "C", "G", "N", "GC"};

static const std::array<std::size_t, Stats::CurveIndexSize> QualityIndex =
    {0, 1, 2, 3, Stats::InvalidIndex, 4, Stats::InvalidIndex};

static const std::array<std::size_t, Stats::CurveIndexSize> ContentIndex = {0, 1, 2, 3, 4, Stats::InvalidIndex, 5};

constexpr std::array<unsigned char, 256> Stats::BaseIndexLookup;
constexpr std::array<signed char, 256> Stats::BaseValueLookup;

// We should use direct initialization syntax
Stats::Stats(Options* opt, bool isRead2, int guessedCycles, int bufferMargin){
    mOptions = opt;
    mIsRead2 = isRead2;
    mReads = 0;
    mLengthSum = 0;

    mEvaluatedSeqLen = mOptions->seqLen1;
    if(mIsRead2)
        mEvaluatedSeqLen = mOptions->seqLen2;

    if(guessedCycles == 0) {
        guessedCycles = mEvaluatedSeqLen;
    }

    mCycles = guessedCycles;
    mBases = 0;
    mQ20Total = 0;
    mQ30Total = 0;
    mQ40Total = 0;
    summarized = false;
    mKmerMin = 0;
    mKmerMax = 0;

    // extend the buffer to make sure it's long enough
    mBufLen = guessedCycles + bufferMargin;

    mQ20Bases.fill(0);
    mQ30Bases.fill(0);
    mBaseContents.fill(0);

    mCycleQ30Bases.resize(Stats::BaseCount, mBufLen);
    mCycleQ20Bases.resize(Stats::BaseCount, mBufLen);
    mCycleBaseContents.resize(Stats::BaseCount, mBufLen);
    mCycleBaseQual.resize(Stats::BaseCount, mBufLen);

    mCycleTotalBase.assign(mBufLen, 0);
    mCycleTotalQual.assign(mBufLen, 0);

    // allocate k-mer counting buffer and build label cache
    mKmer.fill(0);
    for (std::size_t i = 0; i < (Stats::KmerCount >> 4); ++i) {
        auto prefix = kmer3(static_cast<int>(i));
        for (std::size_t j = 0; j < 16; ++j) {
            auto idx = (i << 4) | j;
            mKmerLabels[idx] = prefix + kmer2(static_cast<int>(j));
        }
    }

    initOverRepSeq();
}

void Stats::extendBuffer(int newBufLen) {
    if (newBufLen <= mBufLen) {
        return;
    }

    mCycleQ30Bases.resize(Stats::BaseCount, newBufLen);
    mCycleQ20Bases.resize(Stats::BaseCount, newBufLen);
    mCycleBaseContents.resize(Stats::BaseCount, newBufLen);
    mCycleBaseQual.resize(Stats::BaseCount, newBufLen);

    mCycleTotalBase.resize(newBufLen, 0);
    mCycleTotalQual.resize(newBufLen, 0);

    mBufLen = newBufLen;
}

void Stats::summarize(bool forced) {
    if(summarized && !forced)
        return;

    // first get the cycle and count total bases
    for (std::size_t c = 0; c < mCycleTotalBase.size(); c++) {
        mBases += mCycleTotalBase[c];
        if (mCycleTotalBase[c] == 0) { 
            mCycles = c;
            break;
        }
    }

    if (!mCycleTotalBase.empty() && mCycleTotalBase.back() > 0) {
        mCycles = mBufLen;
    }

    // Q20, Q30, base content
    for (int ci = 0; ci < mCycles; ci++) {
        for (std::size_t i = 0; i < Stats::BaseCount; i++) {
            mQ20Bases[i] += mCycleQ20Bases(ci, i);
            mQ30Bases[i] += mCycleQ30Bases(ci, i);
            mBaseContents[i] += mCycleBaseContents(ci, i);
        }
    }

    mQ20Total = std::accumulate(mQ20Bases.begin(), mQ20Bases.end(), 0L);
    mQ30Total = std::accumulate(mQ30Bases.begin(), mQ30Bases.end(), 0L);

    // This reuses the values from the original loop, not sure what they represent though.
    mQ40Total = std::accumulate(mBaseQualHistogram.begin() + 73,
                                mBaseQualHistogram.begin() + 127, 0L);

    // quality curve for mean qual
    std::vector<double> meanQualCurve(mCycles);
    for(int c=0; c<mCycles; c++) {
        meanQualCurve[c] = (double)mCycleTotalQual[c] / (double)mCycleTotalBase[c];
    }

    // quality curves and base content curves for different nucleotides
    constexpr std::size_t QualityBaseCount = Stats::QualityCurveSize - 1;
    for (std::size_t i = 0; i < QualityBaseCount; i++) {
        const auto& key = QualityOrder[i];
        char base = QualityNames[i][0]; // "A", "T", "C", "G"
        std::size_t baseIdx = baseIndex(base);

        std::vector<double> qualCurve(mCycles);
        std::vector<double> contentCurve(mCycles);

        for (int ci = 0; ci < mCycles; ci++) {
            if (mCycleBaseContents(ci, baseIdx) == 0) {
                qualCurve[ci] = meanQualCurve[ci];
            }
            else {
                qualCurve[ci] = (double)mCycleBaseQual(ci, baseIdx) / (double)mCycleBaseContents(ci, baseIdx);
            }
            contentCurve[ci] = (double)mCycleBaseContents(ci, baseIdx) / (double)mCycleTotalBase[ci];
        }
        mQualityCurves[qualityCurveIndex(key)] = std::move(qualCurve);
        mContentCurves[contentCurveIndex(key)] = std::move(contentCurve);
    }
    mQualityCurves[qualityCurveIndex(CurveKey::Mean)] = std::move(meanQualCurve);

    // GC content curve
    std::vector<double> gcContentCurve(mCycles);
    auto gBase = baseIndex('G');
    auto cBase = baseIndex('C');

    for (int ci = 0; ci < mCycles; ci++) {
        gcContentCurve[ci] = (double)(mCycleBaseContents(ci, gBase) + mCycleBaseContents(ci, cBase)) / (double)mCycleTotalBase[ci];
    }
    mContentCurves[contentCurveIndex(CurveKey::GC)] = std::move(gcContentCurve);

    auto minmax = std::minmax_element(mKmer.begin(), mKmer.end());
    mKmerMin = *minmax.first;
    mKmerMax = *minmax.second;
    mKmerTotal = std::accumulate(mKmer.begin(), mKmer.end(), 0L);

    summarized = true;
}

int Stats::getMeanLength() {
    if(mReads == 0)
        return 0.0;
    else
        return mLengthSum/mReads;
}

void Stats::statRead(Read* r) {
    int len = r->length();

    mLengthSum += len;

    if(mBufLen < len) {
        extendBuffer(std::max(len + 100, (int)(len * 1.5)));
    }
    const char* seqstr = r->mSeq->c_str();
    const char* qualstr = r->mQuality->c_str();

    int kmer = 0;
    bool needFullCompute = true;
    const char q20 = '5';
    const char q30 = '?';

    auto& q20Vec = mCycleQ20Bases.values();
    auto& q30Vec = mCycleQ30Bases.values();
    auto& contentVec = mCycleBaseContents.values();
    auto& qualVec = mCycleBaseQual.values();

    for(int i=0; i<len; i++) {
        char base = seqstr[i];
        char qual = qualstr[i];
        const auto baseIdx = baseIndex(base);

        mBaseQualHistogram[qual]++;

        const int isQ20 = (qual >= q20) ? 1 : 0;
        const int isQ30 = (qual >= q30) ? 1 : 0;
        
        const auto offset = (static_cast<std::size_t>(i) * Stats::BaseCount) + baseIdx;

        q20Vec[offset] += isQ20;
        q30Vec[offset] += isQ30;

        contentVec[offset]++;
        qualVec[offset] += (qual-33);

        mCycleTotalBase[i]++;
        mCycleTotalQual[i] += (qual-33);

        if(base == 'N'){
            needFullCompute = true;
            continue;
        }

        // 5 bases required for kmer computing
        if(i<4)
            continue;

        // calc 5 KMER
        // use Stats::KmerMask to maintain a 10-bit rolling k-mer
        if(!needFullCompute){
            int val = base2val(base);
            if(val < 0){
                needFullCompute = true;
                continue;
            } else {
                kmer = ((kmer<<2) | val) & Stats::KmerMask;
                mKmer[kmer]++;
            }
        } else {
            bool valid = true;
            kmer = 0;
            for(int k=0; k<5; k++) {
                int val = base2val(seqstr[i - 4 + k]);
                if(val < 0) {
                    valid = false;
                    break;
                }
                kmer = ((kmer<<2) | val) & Stats::KmerMask;
            }
            if(!valid) {
                needFullCompute = true;
                continue;
            } else {
                mKmer[kmer]++;
                needFullCompute = false;
            }
        }

    }

    // do overrepresentation analysis for 1 of every 100 reads
    if(mOptions->overRepAnalysis.enabled) {
        if(mReads % mOptions->overRepAnalysis.sampling == 0) {
            const int steps[5] = {10, 20, 40, 100, std::min(150, mEvaluatedSeqLen-2)};
            for(int s=0; s<5; s++) {
                int step = steps[s];
                for(int i=0; i<len-step; i++) {
                    std::string seq = r->mSeq->substr(i, step);
                    if(mOverRepSeq.count(seq)>0) {
                        mOverRepSeq[seq]++;
                        for(int p = i; p < seq.length() + i && p < mEvaluatedSeqLen; p++) {
                            mOverRepSeqDist[seq][p]++;
                        }
                        i+=step;
                    }
                }
            }
        }
    }

    mReads++;
}

auto Stats::base2val(char base) noexcept -> std::int8_t {
    return BaseValueLookup[static_cast<signed char>(base)];
}

auto Stats::baseIndex(char base) noexcept -> std::size_t {
    return BaseIndexLookup[static_cast<unsigned char>(base)];
}

auto Stats::qualityCurveIndex(CurveKey key) noexcept -> std::size_t {
    return QualityIndex[static_cast<std::size_t>(key)];
}

auto Stats::contentCurveIndex(CurveKey key) noexcept -> std::size_t {
    return ContentIndex[static_cast<std::size_t>(key)];
}

int Stats::getCycles() {
    if(!summarized)
        summarize();
    return mCycles;
}

long Stats::getReads() {
    if(!summarized)
        summarize();
    return mReads;
}

long Stats::getBases() {
    if(!summarized)
        summarize();
    return mBases;
}

long Stats::getQ20() {
    if(!summarized)
        summarize();
    return mQ20Total;
}

long Stats::getQ30() {
    if(!summarized)
        summarize();
    return mQ30Total;
}

long Stats::getQ40() {
    if(!summarized)
        summarize();
    return mQ40Total;
}

auto Stats::getQualHist() const -> const std::array<long, 128>& {
    return mBaseQualHistogram;
}

long Stats::getGCNumber() {
    if(!summarized)
        summarize();
    return mBaseContents[baseIndex('G')] + mBaseContents[baseIndex('C')];
}

void Stats::print() {
    if(!summarized) {
        summarize();
    }
    std::cerr << "total reads: " << mReads << "\n";
    std::cerr << "total bases: " << mBases << "\n";
    std::cerr << "Q20 bases: " << mQ20Total << "(" << (mQ20Total*100.0)/mBases << "%)" << "\n";
    std::cerr << "Q30 bases: " << mQ30Total << "(" << (mQ30Total*100.0)/mBases << "%)" << "\n";
    std::cerr << "Q40 bases: " << mQ40Total << "(" << (mQ40Total*100.0)/mBases << "%)" << "\n";
}

void Stats::reportJson(std::ostream& ofs, const std::string& padding) {
    ofs << "{" << "\n";

    ofs << padding << "\t" << "\"total_reads\": " << mReads << "," << "\n";
    ofs << padding << "\t" << "\"total_bases\": " << mBases << "," << "\n";
    ofs << padding << "\t" << "\"q20_bases\": " << mQ20Total << "," << "\n";
    ofs << padding << "\t" << "\"q30_bases\": " << mQ30Total << "," << "\n";
    ofs << padding << "\t" << "\"q40_bases\": " << mQ40Total << "," << "\n";
    ofs << padding << "\t" << "\"total_cycles\": " << mCycles << "," << "\n";

    // quality curves
    ofs << padding << "\t" << "\"quality_curves\": {" << "\n";

    for (const auto key : QualityOrder) {
        const auto curveIndex = qualityCurveIndex(key);
        const auto& curve = mQualityCurves[curveIndex];
        const auto* name = QualityNames[curveIndex]; // A helper for this might be useful

        ofs << padding << "\t\t" << "\"" << name << "\":[";
        const auto curveSize = curve.size();
        for(std::size_t c = 0; c < curveSize; c++) {
            ofs << curve[c];
            // not the end
            if(c - 1 != curveSize)
                ofs << ",";
        }
        ofs << "]";
        // not the end;
        if(key != QualityOrder.back())
            ofs << ",";
        ofs << "\n";
    }
    ofs << padding << "\t" << "}," << "\n";

    // content curves
    ofs << padding << "\t" << "\"content_curves\": {" << "\n";

    for (const auto key : ContentOrder) {
        const auto curveIndex = contentCurveIndex(key);
        const auto& curve = mContentCurves[curveIndex];
        const auto* name = ContentNames[curveIndex];

        ofs << padding << "\t\t" << "\"" << name << "\":[";
        const auto curve_size = curve.size();
        for(std::size_t c = 0; c < curve_size; c++) {
            ofs << curve[c];
            // not the end
            if(c != curve_size - 1)
                ofs << ",";
        }
        ofs << "]";
        // not the end;
        if(key != ContentOrder.back())
            ofs << ",";
        ofs << "\n"; 
    }
    ofs << padding << "\t" << "}," << "\n";

    // KMER counting
    ofs << padding << "\t" << "\"kmer_count\": {" << "\n";
    for(int i=0; i<64; i++) {
        std::string first = kmer3(i);
        for(int j=0; j<16; j++) {
            int target = (i<<4) + j;
            long count = mKmer[target];
            std::string last = kmer2(j);
            ofs << padding << "\t\t\"" << first << last << "\":" << count;
            if(j != 16-1)
                ofs << ",";
        }
        if(i != 64-1)
            ofs << "," << "\n";
        else
            ofs << "\n";
    }
    ofs << padding << "\t" << "}," << "\n";

    // over represented seqs
    std::unordered_map<std::string, long>::iterator iter;
    bool first = true;
    ofs << padding << "\t" << "\"overrepresented_sequences\": {" << "\n";
    for(iter=mOverRepSeq.begin(); iter!=mOverRepSeq.end(); iter++) {
        std::string seq = iter->first;
        long count = iter->second;
        if(!overRepPassed(seq, count))
            continue;
        if(!first) {
            ofs << "," << "\n";
        } else
            first = false;
        ofs << padding << "\t\t\"" << seq <<  "\":" << count;
    }
    ofs << padding << "\t" << "}" << "\n";

    ofs << padding << "}," << "\n";
}

std::string Stats::list2string(double* list, int size) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        ss << list[i];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

std::string Stats::list2string(double* list, int size, long* coords) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        // coords is 1,2,3,...
        long start = 0;
        if(i>0)
            start = coords[i-1];
        long end = coords[i];

        double total = 0.0;
        for(int k=start; k<end; k++)
            total += list[k];

        // get average
        if(end == start)
            ss << "0.0";
        else
            ss << total / (end - start);
        //ss << list[coords[i]-1];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

std::string Stats::list2string(long* list, int size) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        ss << list[i];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

void Stats::reportHtml(std::ostream& ofs, const std::string& filteringType, const std::string& readName) {
    reportHtmlQuality(ofs, filteringType, readName);
    reportHtmlContents(ofs, filteringType, readName);
    reportHtmlKMER(ofs, filteringType, readName);
    if(mOptions->overRepAnalysis.enabled) {
        reportHtmlORA(ofs, filteringType, readName);
    }
}

bool Stats::overRepPassed(std::string& seq, long count) {
    int s = mOptions->overRepAnalysis.sampling;
    switch(seq.length()) {
        case 10:
            return s * count > 500;
        case 20:
            return s * count > 200;
        case 40:
            return s * count > 100;
        case 100:
            return s * count > 50;
        default:
            return s * count > 20;
    }
}

void Stats::reportHtmlORA(std::ostream& ofs, const std::string& filteringType, const std::string& readName) {
    // over represented seqs
    double dBases = mBases;
    std::unordered_map<std::string, long>::iterator iter;
    int displayed = 0;

    // KMER
    std::string subsection = filteringType + ": " + readName + ": overrepresented sequences";
    std::string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div  id='" << divName << "'>\n";
    ofs << "<div class='sub_section_tips'>Sampling rate: 1 / " << mOptions->overRepAnalysis.sampling << "</div>\n";
    ofs << "<table class='summary_table'>\n";
    ofs << "<tr style='font-weight:bold;'><td>overrepresented sequence</td><td>count (% of bases)</td><td>distribution: cycle 1 ~ cycle " << mEvaluatedSeqLen << "</td></tr>"<<"\n";
    int found = 0;
    for(iter = mOverRepSeq.begin(); iter != mOverRepSeq.end(); iter++) {
        std::string seq = iter->first;
        long count = iter->second;
        if(!overRepPassed(seq, count))
            continue;
        found++;
        double percent = (100.0 * count * seq.length() * mOptions->overRepAnalysis.sampling)/dBases;
        ofs << "<tr>";
        ofs << "<td width='400' style='word-break:break-all;font-size:8px;'>" << seq << "</td>";
        ofs << "<td width='200'>" << count << " (" << to_string(percent) <<"%)</td>";
        ofs << "<td width='250'><canvas id='" << divName << "_" << seq << "' width='240' height='20'></td>";
        ofs << "</tr>" << "\n";
    }
    if(found == 0)
        ofs << "<tr><td style='text-align:center' colspan='3'>not found</td></tr>" << "\n";
    ofs << "</table>\n";
    ofs << "</div>\n";

    // output the JS
    ofs << "<script language='javascript'>" << "\n";
    ofs << "var seqlen = " << mEvaluatedSeqLen << ";" << "\n";
    ofs << "var orp_dist = {" << "\n";
    bool first = true;
    for(iter = mOverRepSeq.begin(); iter != mOverRepSeq.end(); iter++) {
        std::string seq = iter->first;
        long count = iter->second;
        if(!overRepPassed(seq, count))
            continue;

        if(!first) {
            ofs << "," << "\n";
        } else
            first = false;
        ofs << "\t\"" << divName << "_" << seq << "\":[";
        for(std::size_t i = 0; i < mOverRepSeqDist[seq].size(); i++) {
            if(i !=0 )
                ofs << ",";
            ofs << mOverRepSeqDist[seq][i];
        }
        ofs << "]";
    }
    ofs << "\n};" << "\n";

    ofs << "for (seq in orp_dist) {"<< "\n";
    ofs << "    var cvs = document.getElementById(seq);"<< "\n";
    ofs << "    var ctx = cvs.getContext('2d'); "<< "\n";
    ofs << "    var data = orp_dist[seq];"<< "\n";
    ofs << "    var w = 240;"<< "\n";
    ofs << "    var h = 20;"<< "\n";
    ofs << "    ctx.fillStyle='#cccccc';"<< "\n";
    ofs << "    ctx.fillRect(0, 0, w, h);"<< "\n";
    ofs << "    ctx.fillStyle='#0000FF';"<< "\n";
    ofs << "    var maxVal = 0;"<< "\n";
    ofs << "    for(d=0; d<seqlen; d++) {"<< "\n";
    ofs << "        if(data[d]>maxVal) maxVal = data[d];"<< "\n";
    ofs << "    }"<< "\n";
    ofs << "    var step = (seqlen-1) /  (w-1);"<< "\n";
    ofs << "    for(x=0; x<w; x++){"<< "\n";
    ofs << "        var target = step * x;"<< "\n";
    ofs << "        var val = data[Math.floor(target)];"<< "\n";
    ofs << "        var y = Math.floor((val / maxVal) * h);"<< "\n";
    ofs << "        ctx.fillRect(x,h-1, 1, -y);"<< "\n";
    ofs << "    }"<< "\n";
    ofs << "}"<< "\n";
    ofs << "</script>"<< "\n";
}

bool Stats::isLongRead() {
    return mCycles > 300;
}

void Stats::reportHtmlKMER(std::ostream& ofs, const std::string& filteringType, const std::string& readName) {

    // KMER
    std::string subsection = filteringType + ": " + readName + ": KMER counting";
    std::string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div  id='" << divName << "'>\n";
    ofs << "<div class='sub_section_tips'>Darker background means larger counts. The count will be shown on mouse over.</div>\n";
    ofs << "<table class='kmer_table' style='width:680px;'>\n";
    ofs << "<tr>";
    ofs << "<td></td>";
    // the heading row
    for(int h=0; h<16; h++) 
        ofs << "<td style='color:#333333'>" << mKmerLabels[h].substr(3) << "</td>";
    ofs << "</tr>\n";
    // content
    for(int i=0; i<64; i++) {
        ofs << "<tr>";

        ofs << "<td style='color:#333333'>" << mKmerLabels[i<<4].substr(0, 3) << "</td>";
        for(int j=0; j<16; j++) {
            ofs << makeKmerTD(i,j) ;
        }
        ofs << "</tr>\n";
    }
    ofs << "</table>\n";
    ofs << "</div>\n";
}

std::string Stats::makeKmerTD(int i, int j) {
    int target = (i<<4) + j;
    long val = mKmer[target];
    // Retrieve k-mer from cache
    const std::string& kmer = mKmerLabels[target];
    double meanBases = (double)(mKmerTotal+1) / static_cast<double>(Stats::KmerCount);
    double prop = val / meanBases;
    double frac = 0.5;
    if(prop > 2.0) 
        frac = (prop-2.0)/20.0 + 0.5;
    else if(prop< 0.5)
        frac = prop;

    frac = std::max(0.01, std::min(1.0, frac));
    int r = (1.0-frac) * 255;
    int g = r;
    int b = r;
    stringstream ss;
    ss << "<td style='background:#"; 
    if(r<16)
        ss << "0";
    ss<<hex<<r;
    if(g<16)
        ss << "0";
    ss<<hex<<g;
    if(b<16)
        ss << "0";
    ss<<hex<<b;
    ss << dec << "' title='"<< kmer << ": " << val << "\n" << prop << " times as mean value'>";
    ss << kmer << "</td>";
    return ss.str();
}

std::string Stats::kmer3(int val) {
    const char bases[4] = {'A', 'T', 'C', 'G'};
    std::string ret(3, ' ');
    ret[0] = bases[(val & 0x30) >> 4];
    ret[1] = bases[(val & 0x0C) >> 2];
    ret[2] = bases[(val & 0x03)];
    return ret;
}

std::string Stats::kmer2(int val) {
    const char bases[4] = {'A', 'T', 'C', 'G'};
    std::string ret(2, ' ');
    ret[0] = bases[(val & 0x0C) >> 2];
    ret[1] = bases[(val & 0x03)];
    return ret;
}

void Stats::reportHtmlQuality(std::ostream& ofs, const std::string& filteringType, const std::string& readName) {
    // quality
    std::string subsection = filteringType + ": " + readName + ": quality";
    std::string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    std::string title;

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";
    
    std::string colors[5] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(20,20,20,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << "\n";
    std::string json_str = "var data=[";

    long *x = new long[mCycles];
    int total = 0;
    if(!isLongRead()) {
        for(int i=0; i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
    } else {
        const int fullSampling = 40;
        for(int i=0; i<fullSampling && i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
        // down sampling if it's too long
        if(mCycles>fullSampling) {
            double pos = fullSampling;
            while(true){
                pos *= 1.05;
                if(pos >= mCycles)
                    break;
                x[total] = (int)pos;
                total++;
            }
            // make sure lsat one is contained
            if(x[total-1] != mCycles){
                x[total] = mCycles;
                total++;
            }
        }
    }
    // four bases
    for (auto key : QualityOrder) {
        const auto curveIndex = qualityCurveIndex(key);
        const auto* base = QualityNames[curveIndex];

        json_str += "{";
        json_str += "x:[" + list2string(x, total) + "],";
        json_str += "y:[" + list2string(mQualityCurves[curveIndex].data(), total, x) + "],";
        json_str += "name: '" + std::string(base) + "',";
        json_str += "mode:'lines',";
        // colors uses the same indexing as QualityNames so we can reuse the curve's index
        json_str += "line:{color:'" + colors[curveIndex] + "', width:1}\n";
        json_str += "},";
    }
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'position'";
    // use log plot if it's too long
    if(isLongRead()) {
        json_str += ",type:'log'";
    }
    json_str += "}, yaxis:{title:'quality'}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << "\n";

    delete[] x;
}

void Stats::reportHtmlContents(std::ostream& ofs, const std::string& filteringType, const std::string& readName) {
    // content
    std::string subsection = filteringType + ": " + readName + ": base contents";
    std::string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    std::string title;

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";
    
    std::string colors[6] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(255, 0, 0, 1.0)", "rgba(20,20,20,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << "\n";
    std::string json_str = "var data=[";

    long *x = new long[mCycles];
    int total = 0;
    if(!isLongRead()) {
        for(int i=0; i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
    } else {
        const int fullSampling = 40;
        for(int i=0; i<fullSampling && i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
        // down sampling if it's too long
        if(mCycles>fullSampling) {
            double pos = fullSampling;
            while(true){
                pos *= 1.05;
                if(pos >= mCycles)
                    break;
                x[total] = (int)pos;
                total++;
            }
            // make sure lsat one is contained
            if(x[total-1] != mCycles){
                x[total] = mCycles;
                total++;
            }
        }
    }
    // four bases
    for (const auto key : ContentOrder) {
        const auto curveIndex = contentCurveIndex(key);
        const auto* baseName = ContentNames[curveIndex];

        long count = 0;
        if (std::strlen(baseName) == 1) {
            const auto bb =  baseIndex(baseName[0]);
            if (bb < Stats::BaseCount) {
                count = mBaseContents[bb];
            }
        } else {
            count = mBaseContents[baseIndex('G')] + mBaseContents[baseIndex('C')];
        }

        std::string percentage = std::to_string(static_cast<double>(count) * 100.0 / mBases);
        if(percentage.length()>5)
            percentage = percentage.substr(0,5);
        auto name = std::string(baseName) + "(" + percentage + "%)";

        json_str += "{";
        json_str += "x:[" + list2string(x, total) + "],";
        json_str += "y:[" + list2string(mContentCurves[curveIndex].data(), total, x) + "],";
        json_str += "name: '" + name + "',";
        json_str += "mode:'lines',";
        json_str += "line:{color:'" + colors[curveIndex] + "', width:1}\n";
        json_str += "},";
    }
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'position'";
    // use log plot if it's too long
    if(isLongRead()) {
        json_str += ",type:'log'";
    }
    json_str += "}, yaxis:{title:'base content ratios'}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << "\n";

    delete[] x;
}

Stats* Stats::merge(std::vector<Stats*>& list) {
    if(list.size() == 0)
        return NULL;

    //get the most long cycles
    int cycles = 0;
    for(int t=0; t<list.size(); t++) {
        list[t]->summarize();
        cycles = std::max(cycles, list[t]->getCycles());
    }

    Stats* s = new Stats(list[0]->mOptions, list[0]->mIsRead2, cycles, 0);

    // init overrepresented seq maps
    std::unordered_map<std::string, long>::iterator iter;

    const auto sizeA = s->mCycleQ30Bases.cycleLength();
    for(int t=0; t<list.size(); t++) {
        int curCycles =  list[t]->getCycles();
        // merge read number
        s->mReads += list[t]->mReads;
        s->mLengthSum += list[t]->mLengthSum;

        const auto sizeB = list[t]->mCycleQ30Bases.cycleLength();
        const auto limit =
            std::min({static_cast<std::size_t>(cycles), sizeA, sizeB, static_cast<std::size_t>(curCycles)});

        // merge per cycle counting for different bases
        for (std::size_t ci = 0; ci < limit; ci++) {
            for (std::size_t bi = 0; bi < Stats::BaseCount; bi++) {
                s->mCycleQ30Bases(ci, bi) += list[t]->mCycleQ30Bases(ci, bi);
                s->mCycleQ20Bases(ci, bi) += list[t]->mCycleQ20Bases(ci, bi);
                s->mCycleBaseContents(ci, bi) += list[t]->mCycleBaseContents(ci, bi);
                s->mCycleBaseQual(ci, bi) += list[t]->mCycleBaseQual(ci, bi);
            }
        }

        // merge per cycle counting for all bases
        for (std::size_t j = 0; j < limit; j++) {
            s->mCycleTotalBase[j] += list[t]->mCycleTotalBase[j];
            s->mCycleTotalQual[j] += list[t]->mCycleTotalQual[j];
        }

        // merge kMer
        std::transform(s->mKmer.begin(), s->mKmer.end(), list[t]->mKmer.begin(),
                       s->mKmer.begin(), std::plus<long>());

        // merge base/read qual histogram
        for(int i=0; i<128; i++) {
            s->mBaseQualHistogram[i] += list[t]->mBaseQualHistogram[i];
        }

        // merge over rep seq
        for(iter = s->mOverRepSeq.begin(); iter != s->mOverRepSeq.end(); iter++) {
            std::string seq = iter->first;
            s->mOverRepSeq[seq] += list[t]->mOverRepSeq[seq];
            if(s->mIsRead2 != list[t]->mIsRead2)
                std::cerr << t <<seq<< ":" << (s->mIsRead2?2:1 ) << "," << (list[t]->mIsRead2?2:1 ) <<"\n";

            std::size_t sizeA = s->mOverRepSeqDist[seq].size();
            std::size_t sizeB = list[t]->mOverRepSeqDist[seq].size();
            for (std::size_t i = 0; i < sizeA && i < sizeB; i++) {
                s->mOverRepSeqDist[seq][i] += list[t]->mOverRepSeqDist[seq][i];
            }
        }
    }

    s->summarize();

    return s;
}

void Stats::initOverRepSeq() {
    std::unordered_map<std::string, long> overRepSeq;
    if(mIsRead2)
        overRepSeq = mOptions->overRepSeqs2;
    else
        overRepSeq = mOptions->overRepSeqs1;

    std::unordered_map<std::string, long>::iterator iter;
    for(iter = overRepSeq.begin(); iter!=overRepSeq.end(); iter++) {
        std::string seq = iter->first;
        mOverRepSeq[seq] = 0;
        mOverRepSeqDist[seq].assign(mEvaluatedSeqLen, 0);
    }
}

auto Stats::test() -> bool {
    Options opt;
    Stats   stats(&opt, false, 10);
    Read    r("@kmer", "ATCGAATCGA", "+", "IIIIIIIIII");
    stats.statRead(&r);
    stats.summarize(true);

    auto idx = [](const std::string& key) -> int {
        unsigned int val = 0;
        for (char ch : key) {
            int base = Stats::base2val(ch);
            if (base < 0) return -1;
            // Shift val left by 2 bits to make space for the next base,
            // mask with 0x3FC (binary 11 1111 1100) to keep only the last 10 bits (rolling window),
            // then add the new base value in the lowest 2 bits.
            val = ((val << 2) & static_cast<unsigned int>(0x3FC)) | static_cast<unsigned int>(base);
        }
        return static_cast<int>(val);
    };

    struct KeyVal {
        const char* key;
        int         value;
    } expected[] = {
        {"ATCGA", 2},
        {"TCGAA", 1},
        {"CGAAT", 1},
        {"GAATC", 1},
        {"AATCG", 1},
    };

    for (const auto& kvPair : expected) {
        auto index = idx(kvPair.key);
        if (index < 0 || stats.mKmer[index] != kvPair.value) {
            return false;
        }
    }

    auto none = idx("AAAAA");
    if (none < 0 || stats.mKmer[none] != 0) {
        return false;
    }

    if (stats.mKmerMax != 2 || stats.mKmerMin != 0) {
        return false;
    }

    if (stats.mKmerTotal != 6) {
        return false;
    }

    return true;
}