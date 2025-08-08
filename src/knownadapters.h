#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "util.h"

namespace adapters {  // adapters

// Forward declaration
class AhoCorasick;

/**
 * @brief Adapter matching result structure
 */
struct AdapterMatch {
    std::string sequence;    // The matched adapter sequence
    std::string name;        // The adapter name/description
    std::size_t position;    // Position in the read where the match starts
    int         mismatches;  // Number of mismatches (for approximate matching)
};

/**
 * @brief Get the map of known adapter sequences
 * @return Reference to the static map of adapter sequences to names/descriptions
 */
auto getKnown() -> const std::unordered_map<std::string, std::string>&;

/**
 * @brief Match a sequence against known adapters using exact prefix matching
 * @param seq The sequence to check
 * @return The matched adapter sequence if found, empty string otherwise
 *
 * @deprecated: use matchKnownAhoCorasick for better performacne
 */
DEPRECATED("use matchKnownAhoCorasick for better performance")
auto matchKnown(const std::string& seq) -> std::string;

/**
 * @brief Get the Aho-Corasick automaton for efficient multi-pattern matching
 * @return Reference to the singleton Aho-Corasick automaton
 *
 * This function initializes the automaton on first call and returns the same instance on on
 * subsequent calls (lazy initialization).
 */
auto getAhoCorasickMatcher() -> AhoCorasick&;

/**
 * @brief Match a sequence against all known adapters using Aho-Corasick
 * @param seq The sequence to check
 * @param startPos Starting position in the sequence (default: 0)
 * @return Vector of all matches found
 *
 * This is significantly more efficient than matchKnown() as it checks
 * all adapters in a single pass with O(n) complexity.
 */
auto matchKnownAhoCorasick(const std::string& seq, std::size_t startPos = 0)
    -> std::vector<AdapterMatch>;

/**
 * @brief Find the first adapter match in a sequence
 * @param seq The sequence to check
 * @param startPos Starting position in the sequence (default: 0)
 * @return Pair of (found, match) - found is true if a match was found
 *
 * This is optimized for cases where only the first match is needed,
 * allowing early termination.
 */
auto findFirstAdapter(const std::string& seq, std::size_t startPos = 0)
    -> std::pair<bool, AdapterMatch>;

/**
 * @brief Match adapters with approximate matching (allowing mismatches)
 * @param seq The sequence to check
 * @param maxMismatches Maximum number of mismatches allowed
 * @param minMatchLength Minimum length of match required
 * @return Vector of matches with their mismatch counts
 *
 * This function extends Aho-Corasick with approximate matching, useful for handling sequencing
 * errors.
 */
auto matchKnownApproximate(const std::string& seq, int maxMismatches = 2, int minMatchLength = 8)
    -> std::vector<AdapterMatch>;

}  // namespace adapters
