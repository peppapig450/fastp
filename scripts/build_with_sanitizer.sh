#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${SCRIPT_DIR}/.."
cd "${PROJECT_ROOT}"

TARGET="fastp"
BUILD_TYPE="$1"

SANITIZER_FLAGS=""
LDFLAGS=""

case "$BUILD_TYPE" in
    asan)
        SANITIZER_FLAGS="-fsanitize=address -fno-omit-frame-pointer -O1 -g"
        LDFLAGS="-fsanitize=address"
        ;;
    msan)
        SANITIZER_FLAGS="-fsanitize=memory -fno-omit-frame-pointer -O1 -g"
        LDFLAGS="-fsanitize=memory"
        ;;
    tsan)
        SANITIZER_FLAGS="-fsanitize=thread -O1 -g"
        LDFLAGS="-fsanitize=thread"
        ;;
    ubsan)
        SANITIZER_FLAGS="-fsanitize=undefined -O1 -g"
        LDFLAGS="-fsanitize=undefined"
        ;;
    valgrind)
        SANITIZER_FLAGS="-O0 -g -fno-omit-frame-pointer"
        ;;
    *)
        echo "Usage: $0 {asan|msan|tsan|ubsan|valgrind}"
        exit 1
        ;;
esac

echo ">>> Cleaning previous builds..."
make clean

echo ">>> Building with $BUILD_TYPE sanitizer..."

cores="$(( $(getconf _NPROCESSORS_ONLN) / 2 ))"
make -j${cores} CXX="clang++" \
     CXXFLAGS="$SANITIZER_FLAGS -std=c++11 -pthread -I./inc" \
     LD_FLAGS="$LDFLAGS -lisal -ldeflate -lpthread" \
     ${TARGET}
