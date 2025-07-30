#include "sequence.h"

#include <cstddef>
#include <iostream>
#include <ostream>
#include <utility>

namespace {                                 // anon

constexpr unsigned char kLowerMask = 0x20;  // mask for lowercasing
constexpr unsigned char kAtMask    = 0x15;  // mask for A/a <-> T/t
constexpr unsigned char kCgMask    = 0x04;  // mask for C/c <-> G/g
}  // namespace

/* DNA base pairing: A <-> T, C <-> G.
 *
 * ASCII character codes for DNA bases have specific patterns:
 * 'A' (0x41) ^ 'T' (0x54) == 0x15
 * 'C' (0x43) ^ 'G' (0x47) == 0x04
 * And this holds for lowercase too:
 * 'a' (0x61) ^ 't' (0x74) == 0x15
 * 'c' (0x63) ^ 'g' (0x67) == 0x04
 * the appropriate bitmask (0x15 for A/T, 0x04 for C/G).
 *
 * So to compute the DNA complement, we lowercase the base
 * (to handle both cases uniformly), and then XOR it with the appropriate
 * bitmask (0x15 for A/T, 0x04 for C/G).
 */
auto Sequence::complement(char base) noexcept -> char {
  // fold to lowercase once
  const auto lower = static_cast<unsigned char>(base) | kLowerMask;

  const auto mask = (lower == 'a' || lower == 't')   ? kAtMask
                    : (lower == 'c' || lower == 'g') ? kCgMask
                                                     : 0;

  // if mask is zero, it wasn't valid so we return 'N'
  return (mask != 0) ? static_cast<char>(static_cast<unsigned char>(base) ^
                                         static_cast<unsigned char>(mask))
                     : 'N';
}

auto Sequence::reverseComplement() const noexcept -> Sequence {
    const std::size_t length = mStr.length();
    std::string       revCompl(length, '\0');

    auto const reverseBase = length - 1;
    for (std::size_t i = 0; i < length; ++i) {
        revCompl[reverseBase - i] = complement(mStr[i]);
    }

    return Sequence(std::move(revCompl));
}

auto operator<<(std::ostream& os, const Sequence& seq) -> std::ostream& { return os << seq.str(); }

auto Sequence::test() noexcept -> bool {
    {
        Sequence s("AAAATTTTCCCCGGGG");
        Sequence rc = ~s;
        if (rc.str() != "CCCCGGGGAAAATTTT") {
            std::cerr << "Failed in reverseComplement() expect CCCCGGGGAAAATTTT: but got " << rc
                      << '\n';
            return false;
        }

        if (s.str() != "AAAATTTTCCCCGGGG") {
            std::cerr << "Failed in reverseComplement() expect AAAATTTTCCCCGGGG: but got " << ~s
                      << '\n';
            return false;
        }
    }
    {
        Sequence s("aAtT");  // mixed-case round trip
        if ((~s).str() != "AaTt") {
            std::cerr << "Test-2 failed, expected 'AaTt': got " << ~s << '\n';
            return false;
        }
    }
    return true;
}
