cmake_minimum_required(VERSION 3.16)

project(evaluator_benchmark
  VERSION 0.1.0
  LANGUAGES CXX)

include(CTest) # defines BUILD_TESTING option

option(ENABLE_IPO "Enable interprocedural optimization (LTO) if supported" ON)
option(ENABLE_ASAN "Enable AddressSanitizer in Debug builds" OFF)

# Sensible default build type if single-config generator (make/ninja)
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

find_package(Threads REQUIRED)

# Google Benchmark (system package)
find_package(benchmark REQUIRED)

# libisal via pkg-config
find_package(PkgConfig REQUIRED)
pkg_check_modules(ISAL REQUIRED IMPORTED_TARGET libisal)
find_package(Python3 COMPONENTS Interpreter REQUIRED)

# Core library made from your ../src files
add_library(evaluator_core STATIC
    ../src/evaluator.cpp
    ../src/options.cpp
    ../src/read.cpp
    ../src/fastqreader.cpp
    ../src/readpool.cpp
    ../src/fastareader.cpp
    ../src/nucleotidetree.cpp
    ../src/sequence.cpp
    ../src/knownadapters.cpp
)

# If headers in ../src are meant to be public, expose that include dir
target_include_directories(evaluator_core
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/../src
)

# Use per-target C++ standard and warnings
target_compile_features(evaluator_core PUBLIC cxx_std_23)

# Warnings by compiler
target_compile_options(evaluator_core
  PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /permissive->
)

# Config-specific optimization/defines (Release already adds NDEBUG)
# Add any extra tuning flags here if you truly need them.

# Enable IPO/LTO if requested and supported
if(ENABLE_IPO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT ipo_supported OUTPUT ipo_msg)
  if(ipo_supported)
    set_property(TARGET evaluator_core PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(STATUS "IPO not supported: ${ipo_msg}")
  endif()
endif()

# -------- Benchmark executable
add_executable(evaluator_benchmark
    evaluator_benchmark.cpp
    adapter_benchmark.cpp
    ../src/ahocorasick.h
)

target_link_libraries(evaluator_benchmark
  PRIVATE
    evaluator_core
    benchmark::benchmark
    benchmark::benchmark_main
    Threads::Threads
    PkgConfig::ISAL
)

# If you have private includes for the benchmark itself:
target_include_directories(evaluator_benchmark PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# AddressSanitizer opt-in for Debug builds
if(ENABLE_ASAN)
  target_compile_options(evaluator_benchmark PRIVATE
    $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-fsanitize=address -fno-omit-frame-pointer>)
  target_link_options(evaluator_benchmark PRIVATE
    $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-fsanitize=address>)
endif()

# Tests (via CTest)
if(BUILD_TESTING)
  add_test(NAME evaluator_benchmark
           COMMAND evaluator_benchmark
                   --benchmark_min_time=2.0s
                   --benchmark_repetitions=1
                   --benchmark_report_aggregates_only=true)
endif()

# Convenience targets
add_custom_target(run_benchmark
  COMMAND evaluator_benchmark --benchmark-format=json --benchmark_out=benchmark_results.json
  DEPENDS evaluator_benchmark
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMENT "Running Google Benchmark and writing results to benchmark_results.json"
)

add_custom_target(check_benchmarks
  COMMAND evaluator_benchmark --benchmark-format=json --benchmark_out=benchmark_results.json
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/compare_benchmarks.py --tolerance=5 ${CMAKE_CURRENT_SOURCE_DIR}/benchmark_baseline.json benchmark_results.json
  DEPENDS evaluator_benchmark
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMENT "Run benchmark and compare against baseline"
)

if(UNIX AND NOT APPLE) # Linux perf
  add_custom_target(profile_benchmark
    COMMAND perf record -g --call-graph=dwarf -- $<TARGET_FILE:evaluator_benchmark> --benchmark_min_time=10s
    DEPENDS evaluator_benchmark
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Profiling benchmark with perf"
  )
endif()

file(GLOB BENCHMARK_HEADERS CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")

if(BENCHMARK_HEADERS)
  add_library(header_headers INTERFACE)
  target_sources(header_headers INTERFACE ${BENCHMARK_HEADERS})
  target_include_directories(header_headers INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include)
endif()

# Generate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
