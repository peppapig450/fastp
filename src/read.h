#pragma once

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "util.h"

class Read {
public:
    /*
     * Legacy aliases - DO *NOT* delete these until all call-sites are migrated.
     * They are raw pointers (not references) so they *can* can be a rebound
     * after a move/copy.
     */
    std::string* mName {nullptr};
    std::string* mSeq {nullptr};
    std::string* mStrand {nullptr};
    std::string* mQuality {nullptr};

    Read(std::string name,
         std::string seq,
         std::string strand  = "+",
         std::string quality = "",
         bool        phred64 = false) noexcept;

    // We implement custom constructors so that existing call-sites can continue
    // using pointers.
    Read(const Read& other);
    Read(Read&& other) noexcept;

    auto operator=(const Read& rhs) -> Read&;
    auto operator=(Read&& rhs) noexcept -> Read&;

    // Default destructor is fine, RAII cleans up for us
    ~Read() = default;

    auto name() const noexcept -> const std::string& { return mName_; }
    auto seq() const noexcept -> const std::string& { return mSeq_; }
    auto strand() const noexcept -> const std::string& { return mStrand_; }
    auto quality() const noexcept -> const std::string& { return mQuality_; }

    void resize(std::size_t new_len);
    void trimFront(std::size_t trim_len);

    auto length() const noexcept -> std::size_t { return mSeq_.size(); }

    void print(std::ostream& os = std::cerr) const noexcept;
    void printFile(std::ostream& file) const noexcept;

    auto toString() const -> std::string;
    auto toStringWithTag(const std::string& tag) const -> std::string;

    void appendToString(std::string& target) const;
    void appendToStringWithTag(std::string& target, const std::string& tag) const;

    // TODO: we keep these around for now as they are called a lot, callsites need to
    // updated.

    DEPRECATED("Use appendToString(std::string&) instead of pointer version")
    void appendToString(std::string* target) const {
        assert(target != nullptr && "Target string pointer is null");
        appendToString(*target);
    }

    DEPRECATED("Use appendToStringWithTag(std::string&, const std::string&) instead")
    void appendToStringWithTag(std::string* target, const std::string* tag) const {
        assert(target != nullptr && "Target string pointer is null");
        assert(tag != nullptr && "Tag string pointer is null");
        appendToStringWithTag(*target, *tag);
    }

    [[deprecated("Use appendToStringWithTag(std::string&, const std::string&) instead")]]
    void appendToStringWithTag(std::string* target, const char* tag) const {
        assert(target != nullptr && "Target string pointer is null");
        appendToStringWithTag(*target, std::string(tag));
    }

    auto firstIndex() const noexcept -> std::string;
    auto lastIndex() const noexcept -> std::string;
    auto lowQualCount(int quality = 20) const noexcept
        -> int;  // LATER: Could this be a smaller type?
    auto reverseComplement() const noexcept -> Read;

    auto fixMGI() -> bool;

    static auto test() noexcept -> bool;

private:
    void convertPhred64to33() noexcept;
    void refreshLegacyPointers() noexcept;

    // Internal storage
    std::string mName_;
    std::string mSeq_;
    std::string mStrand_;
    std::string mQuality_;
};

class ReadPair {
public:
    explicit ReadPair(std::unique_ptr<Read> left, std::unique_ptr<Read> right) noexcept;

    // Rule-of-zero: compiler generated copy / move / destructor are fine

    auto left() noexcept -> Read* { return mLeft.get(); }
    auto left() const noexcept -> const Read* { return mLeft.get(); }

    auto right() noexcept -> Read* { return mRight.get(); }
    auto right() const noexcept -> const Read* { return mRight.get(); }

    // TODO: these should eventually replace left/right
    auto leftRef() noexcept -> Read& { return *mLeft; }
    auto leftRef() const noexcept -> const Read& { return *mLeft; }

    auto rightRef() noexcept -> Read& { return *mRight; }
    auto rightRef() const noexcept -> Read& { return *mRight; }

    // merge a pair, without consideration of seq error caused false INDEL
    auto fastMerge() const -> std::unique_ptr<Read>;

    static auto test() -> bool;

private:
    std::unique_ptr<Read> mLeft;
    std::unique_ptr<Read> mRight;
};

// TODO: refine this to make sure it works drop in, consider using PACK_SIZE to .reserve
// this could cut down on the amount of syncs needed (i think) which improves performance.

/*
 * ReadPack: A hybrid modern container that bridges legacy C-style array interface with
 * modern C++ RAII memory management.
 *
 * This struct maintains both:
 * - Modern: std::vector<std::unique_ptr<Read>> for safe ownership
 * - Legacy: Raw pointer array view (data/count) for backwards compatibility
 *
 * The legacy view is automatically synchronized whenever the container changes.
 */
struct ReadPack {
    /// Primary storage: owns all Read objects with unique_ptr for automatic cleanup
    std::vector<std::unique_ptr<Read>> ownedReads;

    /// Non-owning raw pointer array - DO NOT delete these pointers.
    /// This is automatically updated by sync() to point into the `ownedReads` vector.
    Read** data {nullptr};

    // Number of elements in the data array (same as ownedReads.size())
    int count {0};

    /// Default constructor: creates empty ReadPack
    ReadPack() = default;

    /// Move constructor from vector: takes ownership of reads and syncs legacy view
    explicit ReadPack(std::vector<std::unique_ptr<Read>> reads)
        : ownedReads {std::move(reads)} {
        sync();  // Ensure legacy data/count are updated
    }

    /// Deleted copy operations: unique_ptr makes this non-copyable by design
    /// This prevents accidental copies that would be expensive and error prone.
    ReadPack(const ReadPack&)                    = delete;
    auto operator=(const ReadPack&) -> ReadPack& = delete;

    /// Move operations for efficient ownership transfer
    ReadPack(ReadPack&&) noexcept                    = default;
    auto operator=(ReadPack&&) noexcept -> ReadPack& = default;

    /// unique_ptr automatically handles cleanup for all read objects
    ~ReadPack() = default;

    auto size() const noexcept -> std::size_t { return ownedReads.size(); }
    auto empty() const noexcept -> bool { return ownedReads.empty(); }

    /// Adds a new Read object to the container
    /// Template allows perfect forwarding of or derived types
    /// Automatically syncs legacy view after insertion
    template <typename ReadType>
    void push_back(ReadType&& read) {
        ownedReads.emplace_back(std::forward<ReadType>(read));
        sync();
    }

    /// Access Read object by index (mutable version)
    /// Returns raw pointer for compatibility with legacy code
    auto operator[](std::size_t idx) -> Read* { return ownedReads[idx].get(); }

    /// Access Read object by index (const version)
    auto operator[](std::size_t idx) const -> const Read* { return ownedReads[idx].get(); }

    /// Type aliases for clean definitions
    using iterator       = std::vector<std::unique_ptr<Read>>::iterator;
    using const_iterator = std::vector<std::unique_ptr<Read>>::const_iterator;

    /// Iterator support for loops and STL algorithms
    auto begin() -> iterator { return ownedReads.begin(); }
    auto end() -> iterator { return ownedReads.end(); }
    auto begin() const -> const_iterator { return ownedReads.begin(); }
    auto end() const -> const_iterator { return ownedReads.end(); }

    // NOTE: this may become a performance bottleneck, if that is the case batching syncs
    // is probably a good idea.

    /**
     * Synchronizes the legacy data/count interface with the modern reads vector.
     *
     * This method:
     * 1. Resizes internal raw pointer storage to match ownedReads.size()
     * 2. Extracts raw pointers from each unique_ptr in reads
     * 3. Updates data to point to the raw pointer array
     * 4. Updates count to match the container size
     *
     * Called automatically after any modification (constructor, push_back, etc.)
     *
     * WARNING: After calling, any existing pointers to 'data' may be
     * invalidated if the internal storage reallocates.
     */
    void sync() {
        const auto reads_size = ownedReads.size();

        // Ensure backing storage can hold all reads
        mRawPtrs.resize(reads_size);

        // Extract raw pointers from each unique_ptr
        for (std::size_t i = 0; i < reads_size; ++i) {
            mRawPtrs[i] = ownedReads[i].get();
        }

        // Update legacy interface to point to our raw pointer array
        data  = mRawPtrs.data();
        count = static_cast<int>(reads_size);
    }

private:
    // TODO: in the future this pattern should be migrated away from a container of raw
    // pointers. Storing raw pointers in a vector (or any container) has much poorer cache
    // locality than storing objects directly. Objects are stored in contiguous memory, improving
    // spatial locality, making it easier for the CPU to prefetch.

    /// Backing storage for raw pointers exposed via 'data' field
    /// This vector holds non-owning pointers extracted from 'ownedReads'
    /// Automatically managed via sync - do NOT modify directly
    std::vector<Read*> mRawPtrs;
};