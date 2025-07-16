# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Development Commands

### Building fastp
```bash
# Build regular version
make -j

# Build static version (Linux)
make -j static

# Clean build artifacts
make clean

# Install to system
sudo make install
```

### Development Tools
```bash
# Run unit tests
./fastp test

# Check version
./fastp --version

# Build with sanitizers for debugging
./scripts/build_with_sanitizer.sh {asan|msan|tsan|ubsan|valgrind}
```

### Dependencies
The project depends on external libraries:
- `libisal` - Intel Storage Acceleration Library for compression
- `libdeflate` - Fast deflate compression library
- Standard C++ threading libraries

## Code Architecture

### Core Components

**Main Processing Flow:**
- `main.cpp` - Entry point with command-line parsing
- `options.cpp/h` - Configuration and option management
- `processor.cpp/h` - Main processing coordinator
- `peprocessor.cpp/h` - Paired-end read processing
- `seprocessor.cpp/h` - Single-end read processing

**I/O and Data Structures:**
- `fastqreader.cpp/h` - FASTQ file reading with threading
- `writer.cpp/h` and `writerthread.cpp/h` - Multi-threaded output writing
- `read.cpp/h` and `sequence.cpp/h` - Read and sequence data structures
- `readpool.cpp/h` - Memory pool for read objects

**Processing Modules:**
- `adaptertrimmer.cpp/h` - Adapter detection and trimming
- `filter.cpp/h` - Quality and length filtering
- `basecorrector.cpp/h` - Base correction for overlapping PE reads
- `polyx.cpp/h` - PolyX tail trimming (polyG, polyA, etc.)
- `duplicate.cpp/h` - Duplication detection and removal
- `umiprocessor.cpp/h` - UMI (Unique Molecular Identifier) processing

**Analysis and Reporting:**
- `stats.cpp/h` - Statistics collection and calculation
- `evaluator.cpp/h` - Quality evaluation and metrics
- `htmlreporter.cpp/h` - HTML report generation
- `jsonreporter.cpp/h` - JSON report generation
- `overlapanalysis.cpp/h` - Read overlap analysis for PE data

**Utilities:**
- `matcher.cpp/h` - Pattern matching for adapter detection
- `nucleotidetree.cpp/h` - Tree structure for sequence analysis
- `threadconfig.cpp/h` - Threading configuration
- `util.h` and `common.h` - Utility functions and common definitions

### Key Design Patterns

**Multi-threading Architecture:**
- Uses producer-consumer pattern with `spsc_ring_buffer.h`
- Separate threads for reading, processing, and writing
- Thread-safe statistics collection

**Memory Management:**
- Object pooling for reads to reduce allocation overhead
- Efficient string handling for sequences
- RAII patterns throughout

**Performance Optimizations:**
- The codebase is highly optimized for speed
- Recent commits show performance improvements in stats calculation and sequence matching
- Uses Intel ISA-L and libdeflate for fast compression

## Testing

The project includes a built-in unit test framework:
```bash
./fastp test
```

Test data is available in the `testdata/` directory with sample R1.fq and R2.fq files.

## Key Features Implementation

- **Adapter Detection**: Automatic detection for most common adapters, manual specification supported
- **Quality Filtering**: Multiple quality metrics with configurable thresholds
- **Multi-format Support**: Handles single-end, paired-end, and interleaved FASTQ
- **Reporting**: Comprehensive HTML and JSON reports with visualizations
- **Batch Processing**: `parallel.py` script for processing multiple files

## Working with the Code

- The codebase follows standard C++11 practices
- All headers are in `src/` (no separate `inc/` directory used)
- Uses standard make build system
- Code is performance-critical - be careful with changes that might affect speed
- Recent work has focused on optimizing stats calculation and sequence matching algorithms