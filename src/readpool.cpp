#include "readpool.h"
#include "util.h"
#include <atomic>
#include <memory.h>
#include <unistd.h>
#include "common.h"

ReadPool::ReadPool(Options* opt){
    mOptions = opt;

    // Calculate per-buffer limit
    unsigned long totalLimit = static_cast<unsigned long>(PACK_SIZE) * PACK_IN_MEM_LIMIT;
    if (mOptions->interleavedInput) {
        totalLimit *= 2;
    }

    static constexpr unsigned long kMinPerBufferLimit = 1024UL;
    mPerBufferLimit = max(kMinPerBufferLimit, totalLimit / mOptions->thread);

    initBufferLists();
}

ReadPool::~ReadPool() {
    cleanup();
}

auto ReadPool::input(int tid, Read* data) -> bool {
    // Fast path: try non-blocking produce first
    // This is the most common case and avoids size() calculation
    if (mBufferLists[tid]->tryProduce(data)) {
        // This is for debugging ONLY
        mTotalProduced.fetch_add(1, std::memory_order_relaxed);
        return true;
    }

    // Slow path: buffer might be full, check size limit
    if (mBufferLists[tid]->size() > mPerBufferLimit) {
        return false;   // This thread's buffer is at capacity
    }

    // Buffer has space but might have temporary contention
    // Use blocking produce (spins briefly if needed)
    mBufferLists[tid]->produce(data);
    mTotalProduced.fetch_add(1, std::memory_order_relaxed); // Also for debug 
    return true;
}

void ReadPool::cleanup() {
    // NOTE: This can be switched to using RAII if mBufferLists is a vector
    for (int tid = 0; tid < mOptions->thread; ++tid) {
        auto* list = mBufferLists[tid];
        while (list->canBeConsumed()) {
            Read* r = list->consume();
            mTotalConsumed.fetch_add(1, std::memory_order_relaxed); // Debugging
            delete r;
        }
        delete list;
    }
    delete[] mBufferLists;
    mBufferLists = nullptr;
}

void ReadPool::initBufferLists() {
    mBufferLists = new SingleProducerSingleConsumerList<Read*>*[mOptions->thread];

    for(int tid = 0; tid < mOptions->thread; tid++) {
        static constexpr size_t kMinBufferCapacity = 2048UL;
        size_t bufferCapacity = max(kMinBufferCapacity, mPerBufferLimit);

        mBufferLists[tid] = new SingleProducerSingleConsumerList<Read*>(bufferCapacity);
    }
}

auto ReadPool::size() -> size_t {
    // Lock-free aggregation
    size_t total = 0;
    for (int tid = 0; tid < mOptions->thread; tid++) {
        total += mBufferLists[tid]->size();
    }
    return total;
}

auto ReadPool::getOne() -> Read* {
    // Cache-friendly round-robin biased towards recently active threads
    int startThread = mLastCheckedThread.load(std::memory_order_relaxed);

    // First pass: check recently active threads (likely to have data)
    for (int i = 0; i < min(4, mOptions->thread); i++) {
        int tid = (startThread + i) % mOptions->thread;

        if (mBufferLists[tid]->canBeConsumed()) {
            Read* r = mBufferLists[tid]->consume();
            mLastCheckedThread.store((tid + 1) % mOptions->thread, std::memory_order_relaxed);

            // Debug
            mTotalConsumed.fetch_add(1, std::memory_order_relaxed);
            return r;
        }
    }

    // Second pass: check all remaining threads
    for (int i = 4; i < mOptions->thread; i++) {
        int tid = (startThread + i) % mOptions->thread;

        if (mBufferLists[tid]->canBeConsumed()) {
            Read* r = mBufferLists[tid]->consume();
            mLastCheckedThread.store((tid+ 1) % mOptions->thread, std::memory_order_relaxed);
            mTotalConsumed.fetch_add(1, std::memory_order_relaxed); // for debug only
            return r;
        }
    }

    return nullptr;
}