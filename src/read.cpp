#include "read.h"

#include <algorithm>
#include <cstring>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>

#include "sequence.h"

Read::Read(std::string name,
           std::string seq,
           std::string strand,
           std::string quality,
           bool        phred64) noexcept
    : mName_(std::move(name))
    , mSeq_(std::move(seq))
    , mStrand_(std::move(strand))
    , mQuality_(std::move(quality)) {
    refreshLegacyPointers();
    if (phred64) {
        convertPhred64to33();
    }
}

// Custom copy constructor that copies the underlying storage of one
// Read object to another.
Read::Read(const Read& other)
    : mName_(other.mName_)
    , mSeq_(other.mSeq_)
    , mStrand_(other.mStrand_)
    , mQuality_(other.mQuality_) {
    refreshLegacyPointers();
}

// Custom move constructor that moves the underlying storage values of one
// Read object to another.
Read::Read(Read&& other) noexcept
    : mName_(std::move(other.mName_))
    , mSeq_(std::move(other.mSeq_))
    , mStrand_(std::move(other.mStrand_))
    , mQuality_(std::move(other.mQuality_)) {
    refreshLegacyPointers();
}

// Custom copy on assignment constructor
auto Read::operator=(const Read& rhs) -> Read& {
    // Make sure the object is not being assigned to itself
    if (this != &rhs) {
        mName_    = rhs.mName_;
        mSeq_     = rhs.mSeq_;
        mStrand_  = rhs.mStrand_;
        mQuality_ = rhs.mQuality_;
        refreshLegacyPointers();
    }
    return *this;
}

// Custom move on assignment constructor
auto Read::operator=(Read&& rhs) noexcept -> Read& {
    // Ensure we are not moving the object into itself
    if (this != &rhs) {
        mName_    = std::move(rhs.mName_);
        mSeq_     = std::move(rhs.mSeq_);
        mStrand_  = std::move(rhs.mStrand_);
        mQuality_ = std::move(rhs.mQuality_);
        refreshLegacyPointers();
    }
    return *this;
}

void Read::resize(std::size_t new_len) {
    if (new_len > length()) {
        return;
    }
    mSeq_.resize(new_len);
    mQuality_.resize(new_len);
}

void Read::trimFront(std::size_t trim_len) {
    // Ensure we do not trim more than we have to
    trim_len = std::min(trim_len, length());
    mSeq_.erase(0, trim_len);
    mQuality_.erase(0, trim_len);
}

void Read::print(std::ostream& os) const noexcept {
    os << mName_ << '\n' << mSeq_ << '\n' << mStrand_ << '\n' << mQuality_ << '\n';
}

// This writes the output of print to any ostream which can include a file, cout, etc.
void Read::printFile(std::ostream& file) const noexcept { print(file); }

auto Read::toString() const -> std::string {
    return mName_ + "\n" + mSeq_ + "\n" + mStrand_ + "\n" + mQuality_ + "\n";
}

auto Read::toStringWithTag(const std::string& tag) const -> std::string {
    return mName_ + " " + tag + "\n" + mSeq_ + "\n" + mStrand_ + "\n" + mQuality_ + "\n";
}

// TODO: this can probably be moved to the header
void Read::appendToString(std::string& target) const { target += toString(); }

void Read::appendToStringWithTag(std::string& target, const std::string& tag) const {
    target += toStringWithTag(tag);
}

auto Read::lowQualCount(int quality) const noexcept -> int {
    const auto quality_threshold = static_cast<char>(quality + 33);
    return static_cast<int>(std::count_if(mQuality_.begin(), mQuality_.end(), [&](char qual) {
        return qual < quality_threshold;
    }));
}

auto Read::firstIndex() const noexcept -> std::string {
    const auto        plusPos = mName_.rfind('+');
    const std::size_t endPos  = (plusPos != std::string::npos) ? plusPos : mName_.size();

    const auto colonPos = mName_.rfind(':', endPos);
    if (colonPos == std::string::npos) {
        return "";
    }

    return mName_.substr(colonPos + 1, endPos - colonPos - 1);
}

auto Read::lastIndex() const noexcept -> std::string {
    const auto pos = mName_.find_last_of(":+");
    if (pos == std::string::npos) {
        return {};
    };

    return mName_.substr(pos + 1);
}

auto Read::reverseComplement() const noexcept -> Read {
    Sequence    reverseComplSeq(mSeq_);
    Sequence    reverseCompl = reverseComplSeq.reverseComplement();
    std::string reverseQual(mQuality_.rbegin(), mQuality_.rend());

    return {mName_, reverseCompl.str(), mStrand_, std::move(reverseQual)};
}

auto Read::fixMGI() -> bool {
    const auto nameSize = mName_.size();
    const auto nameBack = mName_.back();

    if (nameSize >= 2 && (nameBack == '1' || nameBack == '2') && mName_[nameSize - 2] == '/') {
        mName_.insert(mName_.end() - 2, ' ');
        refreshLegacyPointers();
        return true;
    }

    return false;
}

void Read::convertPhred64to33() noexcept {
    constexpr int phredOffset = 64 - 33;

    for (auto& c : mQuality_) {
        c = static_cast<char>(std::max(33, c - phredOffset));
    }
}

void Read::refreshLegacyPointers() noexcept {
    mName    = &mName_;
    mSeq     = &mSeq_;
    mStrand  = &mStrand_;
    mQuality = &mQuality_;
}

auto Read::test() noexcept -> bool {
    Read   r("@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA",
           "CTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCA"
             "TTCTATGATGTAGTTTTCCTTAGGAGGACATTTTTTACATGAAATTATTAACCTAAATAGAGTTGATC",
           "+",
           "AAAAA6EEEEEEEEEEEEEEEEE#EEEEEEEEEEEEEEEEE/"
             "EEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<"
             "EEEEAEEEEEEEEEEEEEEEAEEE/EEEEEEEEEEAAEAEAAEEEAEEAA");
    string idx = r.lastIndex();
    return idx == "GGTCCCGA";
}

ReadPair::ReadPair(std::unique_ptr<Read> left, std::unique_ptr<Read> right) noexcept
    : mLeft(std::move(left))
    , mRight(std::move(right)) {}

auto ReadPair::fastMerge() const -> std::unique_ptr<Read> {
    if (mLeft == nullptr || mRight == nullptr) {
        return nullptr;
    }

    auto reverseComplRight = std::unique_ptr<Read>(new Read(mRight->reverseComplement()));

    const auto& leftSeq   = mLeft->seq();
    const auto& leftQual  = mLeft->quality();
    const auto& rightSeq  = reverseComplRight->seq();
    const auto& rightQual = reverseComplRight->quality();

    const auto leftLen  = leftSeq.length();
    const auto rightLen = rightSeq.length();

    if (leftLen == 0 || rightLen == 0) {
        return nullptr;
    }

    static constexpr int  MinOverlap = 30;
    static constexpr char HighQChar  = '?';  // ASCII 63
    static constexpr char LowQChar   = '0';  // ASCII 48

    if (leftLen < static_cast<std::size_t>(MinOverlap)
        || rightLen < static_cast<std::size_t>(MinOverlap)) {
        return nullptr;
    }

    std::size_t bestOverlapLen = 0;
    std::size_t bestOffset     = 0;
    int         diff           = 0;
    int         lowQualDiff    = 0;

    const std::size_t maxOverlapLen = std::min(leftLen, rightLen);

    // Find the first acceptable overlap (smallest that works)
    for (std::size_t overlap_len = MinOverlap; overlap_len <= maxOverlapLen; ++overlap_len) {
        diff = lowQualDiff       = 0;
        const std::size_t offset = leftLen - overlap_len;

        // Take raw pointers to the underlying data for speed in the loop
        const char* leftSeq_ptr  = leftSeq.data() + offset;
        const char* leftQual_ptr = leftQual.data() + offset;

        std::size_t i = 0;
        for (; i < overlap_len; ++i) {
            if (leftSeq_ptr[i] != rightSeq[i]) {
                ++diff;
                if ((leftQual_ptr[i] >= HighQChar && rightQual[i] <= LowQChar)
                    || (leftQual_ptr[i] <= LowQChar && rightQual[i] >= HighQChar)) {
                    ++lowQualDiff;
                }
                if (diff > lowQualDiff || lowQualDiff >= 3) {
                    break;
                }
            }
        }

        if (i == overlap_len) {
            bestOverlapLen = overlap_len;
            bestOffset     = offset;
            break;
        }
    }

    if (bestOverlapLen == 0) {
        return nullptr;
    }

    // Pre-size destination buffers (leftLen + rightLen - overlap)
    const std::size_t totalLen = bestOffset + rightLen;
    std::string       mergedSeq;
    std::string       mergedQual;
    mergedSeq.reserve(totalLen);
    mergedQual.reserve(totalLen);

    // Left prefix (non-overlapping part of left read)
    mergedSeq.append(leftSeq.data(), bestOffset);
    mergedQual.append(leftQual.data(), bestOffset);

    // Append entire right read (will fix overlap chars below)
    mergedSeq.append(rightSeq);
    mergedQual.append(rightQual);

    // Fix the overlapping region
    for (std::size_t i = 0; i < bestOverlapLen; ++i) {
        const std::size_t pos = bestOffset + i;

        if (leftSeq[pos] != rightSeq[i]) {
            // Chose base/quality from the higher-quality side
            if (leftQual[pos] >= HighQChar && rightQual[i] <= LowQChar) {
                mergedSeq[pos]  = leftSeq[pos];
                mergedQual[pos] = leftQual[pos];
            } else {
                mergedSeq[pos]  = rightSeq[i];
                mergedQual[pos] = rightQual[i];
            }
        } else {
            // Same base: sum Phred qualities (ASCII-33 scale)
            mergedQual[pos] = static_cast<char>(leftQual[pos] + rightQual[i] - 33);
        }
    }

    std::ostringstream name;
    // clang-format off
    name << mLeft->name()
         << " merged offset:" << bestOffset
         << " overlap:"       << bestOverlapLen
         << " diff:"          << diff;
    // clang-format on

    return std::unique_ptr<Read>(
        new Read(name.str(), std::move(mergedSeq), "+", std::move(mergedQual)));
}

auto ReadPair::test() -> bool {
    std::unique_ptr<Read> left(
        new Read("@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA",
                 "TTTTTTCTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATT"
                 "CTTATGAGACTCATAGTCATTCTATGATGTAG",
                 "+",
                 "AAAAA6EEEEEEEEEEEEEEEEE#"
                 "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
                 "EEEEEEEE"));
    std::unique_ptr<Read> right(
        new Read("@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA",
                 "AAAAAACTACACCATAGAATGACTATGAGTCTCATAAGAATGCACTCAACTAGTCATCACTCCTGTGTTT"
                 "TCATAAGAAAAAACAGTGTTAGAGTCCAAGAG",
                 "+",
                 "AAAAA6EEEEE/"
                 "EEEEEEEEEEE#"
                 "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
                 "EEEEEEEE"));

    ReadPair pair(std::move(left), std::move(right));
    auto     merged = pair.fastMerge();
    if (merged == nullptr) {
        return false;
    }

    const auto& mergeResult    = merged->seq();
    const auto  expectedResult = std::string(
        "TTTTTTCTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCAT"
         "TCTATGATGTAGTTTTTT");

    return mergeResult == expectedResult;
}
