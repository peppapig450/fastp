#pragma once

#include <vector>
#include <algorithm>
#include <cstring>

class Matcher {
private:
    // Avoid repeated allocations by using thread-local buffers
    static thread_local std::vector<int> s_leftBuffer;
    static thread_local std::vector<int> s_rightBuffer;

    // Set a maximum comparison length to prevent excessive memory usage
    // NOTE: this might introduce regressions in edge cases but I'm not sure
    static constexpr int MAX_CMP_LENGTH = 100000;

    // The core template provides a unified implementation that handles matching
    // with one insertion allowed, parameterized by whether to return the min diff.
    // A specialized version (Small) exists for short strings, enabling faster
    // execution via stack allocation and reduced overhead.
    // This separation avoids branching and heap-usage in performance-critical loops
    // especially since this function runs inside a hot loop.
}