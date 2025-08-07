#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <queue>
#include <string>
#include <utility>
#include <vector>

namespace adapters {

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
    struct Match {
        std::string pattern;   // The matched adapter sequence
        std::string name;      // The adapter description/name
        std::size_t position;  // Position in text where match starts
        std::size_t length;    // The length of the matched pattern
    };

private:
    // Node in the trie structure
    struct TrieNode {
        // TODO: since DNA only ever has 5-6 characters we could go lower to 6 characters,
        // but we keep it at 256 for now.
        static constexpr int kAlphabetSize = 256;  // Support all ASCII characters

        std::array<std::unique_ptr<TrieNode>, kAlphabetSize> children {};
        std::vector<std::pair<std::string, std::string>> output;  // Pattern-name parts at this node

        TrieNode* failure = nullptr;  // Failure link for the AC algorithm

        TrieNode() = default;
    };

    std::unique_ptr<TrieNode> root;
    bool                      built            = false;  // Track if failure links have been built
    std::size_t               minPatternLength = std::numeric_limits<std::size_t>::max();

    // Helper to get or create child node
    // TODO: come back and rename this to something other than `byteKey` maybe `base`?
    auto getOrCreateChild(TrieNode* node, unsigned char byteKey) -> TrieNode* {
        auto& childPtr = node->children[byteKey];

        if (childPtr == nullptr) {
            // This is equivalent to std::make_unique in C++14 or later
            childPtr = std::unique_ptr<TrieNode>(new TrieNode());
        }

        return childPtr.get();
    }

public:
    AhoCorasick()
        : root(std::unique_ptr<TrieNode>(new TrieNode())) {}

    // Default destructor
    ~AhoCorasick() = default;

    // Disable copy but allow move semantics
    AhoCorasick(const AhoCorasick&)                    = delete;
    auto operator=(const AhoCorasick&) -> AhoCorasick& = delete;
    AhoCorasick(AhoCorasick&&)                         = default;
    auto operator=(AhoCorasick&&) -> AhoCorasick&      = default;

    /**
     * @brief Add a pattern to the automaton
     * @param pattern The adapter sequence to match
     * @param name The adapter description/name
     */
    void addPattern(const std::string& pattern, const std::string& name) {
        if (pattern.empty()) {
            return;
        }

        built            = false;  // Invalidate failure links
        minPatternLength = std::min(minPatternLength, pattern.length());

        TrieNode* current = root.get();
        for (unsigned char base : pattern) {
            current = getOrCreateChild(current, base);
        }
        current->output.emplace_back(pattern, name);
    }

    /**
     * @brief Build failure links for the Aho-Corasick automaton
     * Must be called after all patterns are added and before searching
     */
    void build() {
        if (built) {
            return;
        }

        std::queue<TrieNode*> nodeQueue;

        // Initialize: direct children of root have failure link to root
        for (auto& childNode : root->children) {
            if (childNode != nullptr) {
                childNode->failure = root.get();
                nodeQueue.push(childNode.get());
            }
        }

        // BFS to build failure links
        while (!nodeQueue.empty()) {
            TrieNode* current = nodeQueue.front();
            nodeQueue.pop();

            for (std::size_t symbol = 0; symbol < TrieNode::kAlphabetSize; ++symbol) {
                if (current->children[symbol] == nullptr) {
                    continue;
                }

                TrieNode* childNode = current->children[symbol].get();
                nodeQueue.push(childNode);

                // Find failure link by walking up the failure chain
                TrieNode* fallback = current->failure;
                while (fallback != nullptr && fallback->children[symbol] == nullptr) {
                    fallback = fallback->failure;
                }

                if (fallback != nullptr && fallback->children[symbol] != nullptr) {
                    childNode->failure = fallback->children[symbol].get();
                } else {
                    childNode->failure = root.get();
                }

                // Merge output from failure link (for overlapping patterns)
                if (childNode->failure != nullptr && !childNode->failure->output.empty()) {
                    childNode->output.insert(childNode->output.end(),
                                             childNode->failure->output.begin(),
                                             childNode->failure->output.end());
                }
            }
        }

        built = true;
    }

    /**
     * @brief Search for all pattern occurrences in the text
     * @param text The text to search in
     * @return Vector of all matches found
     */
    auto search(const std::string& text) -> std::vector<Match> {
        if (!built) {
            build();
        }

        std::vector<Match> matches;
        TrieNode*          currentNode = root.get();

        for (std::size_t i = 0; i < text.length(); ++i) {
            auto symbol = static_cast<unsigned char>(text[i]);

            //  Follow failure links until we find a match or reach root
            while (currentNode != root.get() && currentNode->children[symbol] == nullptr) {
                currentNode = currentNode->failure;
            }

            if (currentNode->children[symbol] != nullptr) {
                currentNode = currentNode->children[symbol].get();
            }

            // Report all matches at current position
            for (const auto& match : currentNode->output) {
                const auto& pattern = match.first;
                const auto& name    = match.second;

                matches.push_back({pattern, name, i + 1 - pattern.length(), pattern.length()});
            }
        }

        return matches;
    }

    /**
     * @brief Find the first match in the text (early termination)
     * @param text The text to search in
     * @param startPos Starting position in the text
     * @return Optional match if found
     */
    auto findFirst(const std::string& text, std::size_t startPos = 0) -> std::pair<bool, Match> {
        if (!built) {
            build();
        }

        TrieNode* currentNode = root.get();

        for (std::size_t i = startPos; i < text.length(); ++i) {
            auto symbol = static_cast<unsigned char>(text[i]);

            while (currentNode != root.get() && currentNode->children[symbol] == nullptr) {
                currentNode = currentNode->failure;
            }

            if (currentNode->children[symbol] != nullptr) {
                currentNode = currentNode->children[symbol].get();
            }

            if (!currentNode->output.empty()) {
                const auto& pattern = currentNode->output[0].first;
                const auto& name    = currentNode->output[0].second;
                return {true,
                        {/*pattern=*/pattern,
                         /*name=*/name,
                         /*position=*/i + 1 - pattern.length(),
                         /*length=*/pattern.length()}};
            }
        }

        return {false, {}};
    }

    /**
     * @brief Check if automaton contains any patterns
     */
    auto empty() const noexcept -> bool {
        using ChildPtr = std::unique_ptr<TrieNode>;

        // Check if root has any children
        return std::all_of(root->children.begin(),
                           root->children.end(),
                           [](const ChildPtr& child) -> bool { return child == nullptr; });
    }

    /**
     * @brief Get minimum pattern length (useful for optimization)
     */
    auto getMinPatternLength() const noexcept -> std::size_t { return minPatternLength; }
};

}  // namespace adapters