#include "sequence.h"

#include <cstddef>
#include <iostream>
#include <ostream>
#include <utility>

// TODO: profile this and determine if the cases are fine, if not
// we can switch to using xor with 0x15 for A(0x41) <-> T(0x54),
// and C(0x43) <-> G(0x47) with 0x04.
auto Sequence::complement(char base) noexcept -> char {
    switch (base) {
        case 'A': return 'T'; case 'a': return 't';
        case 'T': return 'A'; case 't': return 'a';
        case 'C': return 'G'; case 'c': return 'g';
        case 'G': return 'C'; case 'g': return 'c';
        default:  return 'N';
    }
}

auto Sequence::reverseComplement() const noexcept -> Sequence  {
    const std::size_t length = mStr.length();
    std::string revCompl(length, '\0');

    auto const reverseBase = length - 1;
    for (std::size_t i = 0; i < length; ++i) {
        revCompl[reverseBase - i] = complement(mStr[i]);
    }

    return Sequence(std::move(revCompl));
}

auto operator<<(std::ostream& os, const Sequence& seq) -> std::ostream& {
    return os << seq.str();
}


auto Sequence::test() noexcept -> bool {
    {
        Sequence s("AAAATTTTCCCCGGGG");
        Sequence rc = ~s;
        if (rc.str() != "CCCCGGGGAAAATTTT") {
            std::cerr << "Failed in reverseComplement() expect CCCCGGGGAAAATTTT: but got " << rc << '\n';
            return false;
        }

        if (s.str() != "AAAATTTTCCCCGGGG") {
            std::cerr << "Failed in reverseComplement() expect AAAATTTTCCCCGGGG: but got " << ~s << '\n';
            return false;
        }
    }
    {
        Sequence s("aAtT"); // mixed-case round trip
        if ((~s).str() != "AaTt") {
            std::cerr << "Test-2 failed, expected 'AaTt': got " << ~s << '\n';
            return false;
        }
    }
    return true;
}