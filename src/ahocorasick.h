#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <queue>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "base/string_ref.hpp"

namespace adapters {

using fastp::base::StringRef;

/**
 * @brief Aho-Corasick automaton for efficient multi-pattern string matching
 *
 * This implementation provides O(n + m + z) time complexity where:
 * - n is the total length of the text being searched
 * - m is the total length of all patterns
 * - z is the number of pattern occurrences
 */
class AhoCorasick final {
public:
    // Lightweight handle to an interned pattern
    using PatternId = std::uint32_t;

    struct PatternEntry {
        StringRef sequence;  // view into getKnown() strings
        StringRef label;     // view into getKnown() strings
    };

    struct Match {
        PatternId   id;        // id of interned (pattern, name)
        std::size_t position;  // Position in text where the match starts
        std::size_t length;    // The length of the matched patterns
    };

private:
    // DNA Alphabet mapper
    enum class Alpha : unsigned char { A = 0, C = 1, G = 2, T = 3, N = 4, Other = 5, Count = 6 };

    using Index                               = std::size_t;
    using NodeIndex                           = std::int32_t;
    static constexpr auto        InvalidNode  = static_cast<NodeIndex>(-1);
    static constexpr auto        AlphabetSize = static_cast<Index>(Alpha::Count);
    static constexpr std::size_t AlphaLutSize = 26;
    using AlphaArray                          = std::array<Alpha, AlphaLutSize>;

    // Build a 26 entry lookup table once.
    static auto alphaLut26() -> const AlphaArray& {
        static const AlphaArray lookupTable = {{/* a */ Alpha::A,
                                                /* b */ Alpha::Other,
                                                /* c */ Alpha::C,
                                                /* d */ Alpha::Other,
                                                /* e */ Alpha::Other,
                                                /* f */ Alpha::Other,
                                                /* g */ Alpha::G,
                                                /* h */ Alpha::Other,
                                                /* i */ Alpha::Other,
                                                /* j */ Alpha::Other,
                                                /* k */ Alpha::Other,
                                                /* l */ Alpha::Other,
                                                /* m */ Alpha::Other,
                                                /* n */ Alpha::N,
                                                /* o */ Alpha::Other,
                                                /* p */ Alpha::Other,
                                                /* q */ Alpha::Other,
                                                /* r */ Alpha::Other,
                                                /* s */ Alpha::Other,
                                                /* t */ Alpha::T,
                                                /* u */ Alpha::T,  // treat U as T
                                                /* v */ Alpha::Other,
                                                /* w */ Alpha::Other,
                                                /* x */ Alpha::Other,
                                                /* y */ Alpha::Other,
                                                /* z */ Alpha::Other}};

        return lookupTable;
    }

    static auto getAlphaIndex(unsigned char nucleotide) noexcept -> Alpha {
        constexpr unsigned char asciiLowerMask = 0x20;

        // ASCII lowercase only if A..Z
        // NOTE: this can be made branchless using bit masking but that would be a bit of a
        // premature optimization for now
        const unsigned char lowerCase =
            (nucleotide >= 'A' && nucleotide <= 'Z')
                ? static_cast<unsigned char>(nucleotide | asciiLowerMask)
                : nucleotide;

        // If 'a'..'z', use the 26-entry lookup table, otherwise use Other
        const auto offset = static_cast<unsigned char>(lowerCase - 'a');
        if (offset < static_cast<unsigned char>(AlphaLutSize)) {
            return alphaLut26()[offset];
        }

        return Alpha::Other;
    }

    static auto getAlphaIndex(char nucleotide) noexcept -> Alpha {
        return getAlphaIndex(static_cast<unsigned char>(nucleotide));
    }

    struct AutomatonNode {
        // next state for each alphabet symbol (-1 if not set before build())
        std::array<std::int32_t, AlphabetSize> nextState;

        std::int32_t failureLink {0};  // fallback link
        std::int32_t outputHead {-1};  // index into outputs_ list, -1 if none

        AutomatonNode()
            : nextState() {
            nextState.fill(-1);
        }
    };

    struct OutputLink {
        PatternId    patternId;
        std::int32_t nextLink;  // next in linked list, -1 if none
    };

    std::vector<AutomatonNode> nodes_;     // automaton states (0 = root)
    std::vector<OutputLink>    outputs_;   // all outputs are packed here
    std::vector<PatternEntry>  patterns_;  // registered pattern views
    bool                       isBuilt_ {};
    std::size_t                shortestPatternLength_;

    static auto toIndex(NodeIndex index) -> Index {
        assert(index >= 0);
        return static_cast<Index>(index);
    }
    static auto toIndex(Alpha nucleotide) -> Index { return static_cast<Index>(nucleotide); }

    auto node(NodeIndex index) -> AutomatonNode& { return nodes_[toIndex(index)]; }
    auto node(NodeIndex index) const -> const AutomatonNode& { return nodes_[toIndex(index)]; }

    auto nextNode(NodeIndex state, Index symbolIndex) -> NodeIndex& {
        return nodes_[toIndex(state)].nextState[symbolIndex];
    }

    // Create and return new node index
    auto createNewNode() -> int {
        nodes_.emplace_back();
        return static_cast<int>(nodes_.size() - 1);
    }

    // Add output link from node to pattern
    void appendOutput(NodeIndex nodeIndex, PatternId patternId) {
        auto& node = nodes_[toIndex(nodeIndex)];

        outputs_.push_back({patternId, node.outputHead});
        node.outputHead = static_cast<std::int32_t>(outputs_.size() - 1);
    }

public:
    AhoCorasick()
        : shortestPatternLength_(std::numeric_limits<std::size_t>::max()) {
        nodes_.reserve(1);
        createNewNode();  // root at state 0
    }

    // Default destructor
    ~AhoCorasick() = default;

    // Disable copy but allow move semantics
    AhoCorasick(const AhoCorasick&)                    = delete;
    auto operator=(const AhoCorasick&) -> AhoCorasick& = delete;
    AhoCorasick(AhoCorasick&&)                         = default;
    auto operator=(AhoCorasick&&) -> AhoCorasick&      = default;

    void reserve(std::size_t estimatedNodes, std::size_t estimatedOutputs) {
        if (estimatedNodes > nodes_.size()) {
            nodes_.reserve(estimatedNodes);
        }
        if (estimatedOutputs > outputs_.size()) {
            outputs_.reserve(estimatedOutputs);
        }
        // TODO: since the number of adapters is known we can probably do this better
        patterns_.reserve(estimatedOutputs);
    }

    /**
     * @brief Access an interned pattern entry by id (no copies).
     */
    auto getPatternEntry(PatternId patternId) const -> const PatternEntry& {
        return patterns_[patternId];
    }

    auto getPatterns() const -> const std::vector<PatternEntry>& { return patterns_; }

    /**
     * @brief Add a pattern to the automaton
     * @param pattern The adapter sequence to match
     * @param name The adapter description/name
     */
    void addPattern(const std::string& pattern, const std::string& name) {
        if (pattern.empty()) {
            return;
        }

        isBuilt_               = false;  // Invalidate failure links
        shortestPatternLength_ = std::min(shortestPatternLength_, pattern.length());

        NodeIndex currentNode = 0;  // root
        for (auto nucleotide : pattern) {
            const auto alphabetIndex = getAlphaIndex(static_cast<unsigned char>(nucleotide));
            const auto symbolIndex   = toIndex(alphabetIndex);

            // read onc
            NodeIndex& childNode = nextNode(currentNode, symbolIndex);
            if (childNode == InvalidNode) {
                childNode = createNewNode();
            }

            currentNode = childNode;
        }

        auto newPatternId = static_cast<PatternId>(patterns_.size());
        patterns_.push_back({StringRef(pattern), StringRef(name)});
        appendOutput(currentNode, newPatternId);
    }

    // TODO: this should probably be a orchestrator for two smaller functions since this is
    // essentially doing two different operations, initaliasing the root transitions, and BFS to
    // buid the links
    /**
     * @brief Build failure links for the Aho-Corasick automaton
     * Must be called after all patterns are added and before searching
     */
    void build() {
        if (isBuilt_) {
            return;
        }

        std::queue<int> nodeQueue;

        AutomatonNode& rootNode = nodes_.front();

        // Initialize root transitions and enqueue children
        for (Index alphabetIndex = 0; alphabetIndex < AlphabetSize; ++alphabetIndex) {
            auto childIndex = static_cast<NodeIndex>(rootNode.nextState[alphabetIndex]);
            if (childIndex != InvalidNode) {
                node(childIndex).failureLink = 0;
                nodeQueue.push(childIndex);
            } else {
                rootNode.nextState[alphabetIndex] = 0;  // missing root transitions -> root
            }
        }

        // BFS to build failure links
        while (!nodeQueue.empty()) {
            NodeIndex currentIndex = nodeQueue.front();
            nodeQueue.pop();
            auto& curentNode = node(currentIndex);

            for (Index alphabetIndex = 0; alphabetIndex < AlphabetSize; ++alphabetIndex) {
                NodeIndex childIndex = curentNode.nextState[alphabetIndex];
                if (childIndex != InvalidNode) {
                    nodeQueue.push(childIndex);

                    auto&       childNode   = node(childIndex);
                    const auto& failureNode = node(curentNode.failureLink);

                    childNode.failureLink = failureNode.nextState[alphabetIndex];

                    NodeIndex inheritedNode = node(childNode.failureLink).outputHead;
                    if (inheritedNode != InvalidNode) {
                        Index tail = toIndex(inheritedNode);
                        while (outputs_[tail].nextLink != InvalidNode) {
                            tail = toIndex(outputs_[tail].nextLink);
                        }

                        outputs_[tail].nextLink = childNode.outputHead;
                        childNode.outputHead    = inheritedNode;
                    }
                } else {
                    curentNode.nextState[alphabetIndex] =
                        node(curentNode.failureLink).nextState[alphabetIndex];
                }
            }
        }

        isBuilt_ = true;
    }

    /**
     * @brief Search for all pattern occurrences into a caller-provided vector (reusable
     * buffer).
     * @param text The text to search in
     * @param matches Output vector to be filled; capacity can be reused across calls
     */
    void searchInto(const std::string& text, std::vector<Match>& matches) {
        if (!isBuilt_) {
            build();
        }

        matches.clear();
        matches.reserve(8);

        NodeIndex currentState = 0;
        for (Index i = 0; i < text.size(); ++i) {
            currentState = node(currentState).nextState[toIndex((getAlphaIndex(text[i])))];

            NodeIndex outputIndex = node(currentState).outputHead;
            while (outputIndex != InvalidNode) {
                PatternId           patternId = outputs_[toIndex(outputIndex)].patternId;
                const PatternEntry& entry     = patterns_[patternId];

                matches.push_back({patternId,
                                   i + 1 - entry.sequence.size(),  // record the start of the match
                                   entry.sequence.size()});

                outputIndex = outputs_[toIndex(outputIndex)].nextLink;
            }
        }
    }

    /**
     * @brief Search for all pattern occurrences in the text (convenience wrapper)
     *         Prefer searchInto() when you can reuse a buffer.
     * @param text The text to search in
     * @return Vector of all matches found
     */
    auto search(const std::string& text) -> std::vector<Match> {
        std::vector<Match> matches;
        searchInto(text, matches);
        return matches;
    }

    /**
     * Iterate matches with a callback to avoid building any vector.
     * @param text The text to search in
     * @param onMatch Callable with signature
     *        void(PatternId id, const PatternEntry& entry, std::size_t pos, std::size_t len)
     */
    template <typename MatchCallbackFunction>
    void forEachMatch(const std::string& text, MatchCallbackFunction&& onMatch) {
        using result_t = typename std::result_of<
            MatchCallbackFunction(PatternId, const PatternEntry&, std::size_t, std::size_t)>::type;
        static_assert(std::is_same<void, result_t>::value,
                      "onMatch must be callable as: void(PatternId, const PatternEntry&, "
                      "std::size_t, std::size_t)");

        if (!isBuilt_) {
            build();
        }

        NodeIndex currentState = 0;
        for (Index i = 0; i < text.size(); ++i) {
            currentState = node(currentState).nextState[toIndex((getAlphaIndex(text[i])))];

            NodeIndex outputIndex = node(currentState).outputHead;
            while (outputIndex != InvalidNode) {
                const PatternId     patternId     = outputs_[toIndex(outputIndex)].patternId;
                const PatternEntry& entry         = patterns_[patternId];
                const Index         patternLength = entry.sequence.size();
                const Index         matchPosition = i + 1 - patternLength;

                onMatch(patternId, entry, matchPosition, patternLength);
                outputIndex = outputs_[toIndex(outputIndex)].nextLink;
            }
        }
    }

    /**
     * @brief Find the first match in the text (early termination)
     * @param text The text to search in
     * @param startPos Starting position in the text
     * @return (found, Match) pair; when found==false, Match contents are unspecified
     */
    auto findFirst(const std::string& text, std::size_t startPos = 0) -> std::pair<bool, Match> {
        if (!isBuilt_) {
            build();
        }

        const Index textLength = text.size();
        if (startPos >= textLength) {
            return {false, Match {}};
        }

        NodeIndex currentState = 0;

        // Advance automaton up to startPos
        for (Index pos = 0; pos < startPos; ++pos) {
            const Index symbolIndex = toIndex(getAlphaIndex(text[pos]));
            currentState            = node(currentState).nextState[symbolIndex];
        }

        for (Index pos = startPos; pos < textLength; ++pos) {
            const Index symbolIndex = toIndex(getAlphaIndex(text[pos]));
            currentState            = node(currentState).nextState[symbolIndex];

            NodeIndex outputNodeIndex = node(currentState).outputHead;
            if (outputNodeIndex != InvalidNode) {
                const PatternId patternId =
                    outputs_[toIndex(outputNodeIndex)].patternId;  // first match is enough
                const PatternEntry& patternInfo = patterns_[patternId];
                const Index         matchLength = patternInfo.sequence.size();
                const Index         matchStart  = pos + 1 - matchLength;

                return {true,
                        Match {
                            patternId,    // matched pattern
                            matchStart,   // start position in text
                            matchLength,  // match length
                        }};
            }
        }

        return {false, {}};
    }

    /**
     * @brief Check if automaton contains any patterns
     */
    auto empty() const noexcept -> bool {
        // root has any outgoing edge?
        return patterns_.empty();
    }

    /**
     * @brief Get minimum pattern length (useful for optimization)
     */
    auto getMinPatternLength() const noexcept -> std::size_t { return shortestPatternLength_; }
};

}  // namespace adapters
