// in-memory pooled Reads to reduce new/delete operations
// for each thread, one SISC list is used

#ifndef READ_POOL_H
#define READ_POOL_H

#include <cstdio>
#include <cstdlib>
#include <atomic>
#include "read.h"
#include <memory>
#include <vector>
#include "options.h"
#include "spsc_ring_buffer.h"

using namespace std;

class ReadPool{
public:
    ReadPool(Options* opt);
    ~ReadPool();
    auto input(int tid, Read* r) -> bool;
    auto getOne() -> Read*;
    void initBufferLists();
    void cleanup();
    auto size() -> size_t;

    // Temporary monitoring events
    auto getTotalPressureEvents() const -> size_t;
    void printPressureStats() const;
    void resetPressureStats();

private:
    Options* mOptions;

    std::vector<std::unique_ptr<SingleProducerSingleConsumerList<Read*>>> mBufferLists;
   
    unsigned long mPerBufferLimit;      // Per-buffer limit instead of global

    // These are for debugging ONLY
    mutable atomic<size_t> mTotalProduced{0};
    mutable atomic<size_t> mTotalConsumed{0};

    // Round-robin state for getOne()
    mutable atomic<int> mLastCheckedThread{0}; // Fair round-robin
};

#endif