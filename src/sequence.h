#pragma once

#include <cstddef>
#include <ostream>
#include <string>

/**
 * @brief Lightweight DNA (or RNA-like) sequence wrapper.
 *
 *  * Rule-of-zero: only the string-taking constructor is user-defined;
 *    everything else is compiler-generated.
 *  * Single-pass reverse-complement algorithm.
 *  * Case is preserved ( 'A'→'T', 'a'→'t' ).  Unknown bases map to 'N'.
 */
class Sequence {
public:
    Sequence() = default;
    explicit Sequence(std::string seq) noexcept
        : mStr(std::move(seq)) {}

    // observers
    auto length() const noexcept -> std::size_t { return mStr.length(); }
    auto str() const noexcept -> const std::string& { return mStr; }

    // transformations
    auto reverseComplement() const noexcept -> Sequence;
    auto operator~() const noexcept -> Sequence { return reverseComplement(); }

    static auto test() noexcept -> bool;

private:
    static auto complement(char base) noexcept -> char;

    std::string mStr;
};

// For usage like std::cout << seq;
inline auto operator<<(std::ostream& os, const Sequence& seq) -> std::ostream&;
