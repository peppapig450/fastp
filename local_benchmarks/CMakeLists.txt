cmake_minimum_required(VERSION 3.16)
project(evaluator_benchmark CXX)

# Set language standard & build options
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add your optimization and define flags
add_compile_options(-O2 -DNDEBUG -pthread)

# Include directory for your project headers
include_directories(${CMAKE_SOURCE_DIR}/../src)

# Source files from ../src
set(FASTP_SRCS
    ../src/ahocorasick.h
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

# Benchmark source(s)
set(BENCHMARK_SRCS
    evaluator_benchmark.cpp
)

# Combine all sources
set(SRCS
    ${FASTP_SRCS}
    ${BENCHMARK_SRCS}
)

# Find Google Benchmark
find_package(benchmark REQUIRED)

# Find libisal
find_package(PkgConfig REQUIRED)
pkg_check_modules(ISAL REQUIRED libisal)

include_directories(${ISAL_INCLUDE_DIRS})
link_directories(${ISAL_LIBRARY_DIRS})

# Create executable
add_executable(${PROJECT_NAME} ${SRCS})

# Link Google Benchmark and pthread
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        benchmark
        benchmark_main
        pthread
        ${ISAL_LIBRARIES}
)

# Add custom targets for benchmark & profile (optional)
add_custom_target(run_benchmark
    COMMAND ${PROJECT_NAME} --benchmark-format=json --benchmark_out=benchmark_results.json
    DEPENDS ${PROJECT_NAME}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running Google Benchmark and writing results to benchmark_results.json"
)

add_custom_target(profile_benchmark
    COMMAND perf record -g --call-graph=dwarf -- ${PROJECT_NAME} --benchmark_min_time=10s
    DEPENDS ${PROJECT_NAME}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Profiling benchmark with perf"
)
