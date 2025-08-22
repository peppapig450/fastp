#pragma once

// We are about to do some wizardry with macros and compiler quirks to squeeze performance out of
// the CPU. If you're reading this, proceed with caution before dropping this into arbitrary code.

// Ensure __has_cpp_attribute is safe as not all compilers have it (mostly older ones).
// For example: GCC < 5, Clang < 3.3, MSVC > 2017, version 15.3. If you're targeting those,
// you deserve a medal.
#ifndef __has_cpp_attribute
    #define __has_cpp_attribute(x) 0
#endif

// --- Branch Prediction Macros ---
// This is necessary when processing FASTQ because the data is like kryptonite for branch
// predictors: highly irregular data, random distributions, and lots of character-level checks.
// Compilers are unable to always guess correctly, which leads to dismal performance as branch
// mispredictions cause the entire pipeline to flush which becomes very costly. So we use manual
// hints to help the compiler out.
#if __has_cpp_attribute(likely) && __cplusplus >= 202002L
// Prefer standard C++20 hints when available
    #define LIKELY_BRANCH(x) (x)
    #define UNLIKELY_BRANCH(x) (x)
    #define LIKELY [[likely]]
    #define UNLIKELY [[unlikely]]
#else
    // Fallback, no attributes so we use compiler intrinsic hints if available
    #define LIKELY
    #define UNLIKELY

    #if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) \
        || defined(__INTEL_LLVM_COMPILER)
        #define LIKELY_BRANCH(x) __builtin_expect(!!(x), 1)
        #define UNLIKELY_BRANCH(x) __builtin_expect(!!(x), 0)
    #else
        #define LIKELY_BRANCH(x) (x)
        #define UNLIKELY_BRANCH(x) (x)
    #endif
#endif

// --- Loop unrolling macros ---
// Reminder: loop unrolling ONLY helps in specific cases
// NOTE: MSVC doesn't support a numeric unroll count via __pragma(loop(...))
// so we ignore it
#ifndef UNROLL_LOOP
    #if defined(__clang__) || defined(__INTEL_LLVM_COMPILER)
        #include "preprocessor.hpp"
        #define UNROLL_LOOP(N) _Pragma(TOSTRING(clang loop unroll_count(N)))
    #elif defined(__INTEL_COMPILER)
        #include "preprocessor.hpp"
        #define UNROLL_LOOP(N) _Pragma(TOSTRING(unroll(N)))
    #elif defined(__GNUC__)
        #include "preprocessor.hpp"
        #define UNROLL_LOOP(N) _Pragma(TOSTRING(GCC unroll N))
    #elif defined(_MSC_VER)
        #include "preprocessor.hpp"
        #define UNROLL_LOOP(N) __pragma(loop(unroll))
    #else
        #define UNROLL_LOOP(N)
    #endif
#endif

// Macro to apply compiler-specific restrict qualifiers for pointers.
//
// `RESTRICT` hints to the compiler that pointers annotated with it do NOT alias with other
// pointers, allowing for more aggressive optimizations (e.g., vectorization, better instruction
// scheduling).
//
// WARNING: Only apply `RESTRICT` to pointers you can GUARANTEE do NOT alias with other pointers.
// Misusing restrict leads to undefined behavior and quiet miscompilations.
//
// TLDR: Use this on local, non-overlapping buffers (e.g., scratch space) for performance wins.
// If you’re unsure, do not use `RESTRICT`. Bad restrict is worse than no restrict.
#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) \
    || defined(__INTEL_LLVM_COMPILER)
    #define RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define RESTRICT __restrict
#else
    #define RESTRICT
#endif
