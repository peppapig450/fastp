#include "readpool.h"
#include "spsc_ring_buffer.h"
#include "util.h"
#include <atomic>
#include <memory.h>
#include <memory>
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

    // Leverage relaxed capacity checks due to improved backoff efficiency.
    // Inputs are only rejected when usage exceeds 2 the per-buffer limit,
    // allowing higher tolerance before dropping data.
    if (mBufferLists[tid]->size() > mPerBufferLimit * 2) {
        return false;   // Drop only when severely over capacity
    }

    // Buffer has space but might have temporary contention
    // Use blocking produce (spins briefly if needed)
    mBufferLists[tid]->produce(data);
    mTotalProduced.fetch_add(1, std::memory_order_relaxed); // Also for debug 
    return true;
}

void ReadPool::cleanup() {
    for (auto& buffer : mBufferLists) {
        while (buffer->canBeConsumed()) {
            Read *r = buffer->consume();
            mTotalConsumed.fetch_add(1, std::memory_order_relaxed);
            delete r;
        }
    }
    mBufferLists.clear();
}

void ReadPool::initBufferLists() {
    // More efficient backoff means we can use larger buffers for throughput
    mBufferLists.reserve(mOptions->thread);

    for (int tid = 0; tid < mOptions->thread; tid++) {
        // TODO: figure out whether the 8MB default is really needed
        static constexpr size_t kMinBufferCapacity = 4096UL;
        size_t bufferCapacity = std::max(kMinBufferCapacity, mPerBufferLimit * 2);

        // use SPSCRingBuffer use its 8MiB default if our calculation is smaller
        if (bufferCapacity < 1024 * 1024) { 
            bufferCapacity = 0; // Use default 8MiB from SPSCRingBuffer
        }

        mBufferLists.emplace_back(
            // Ensure unique for perfect forwarding
            std::unique_ptr<SingleProducerSingleConsumerList<Read*>>(
                new SingleProducerSingleConsumerList<Read*>(bufferCapacity)
            )
        );
    }
}

// Monitoring methods
auto ReadPool::getTotalPressureEvents() const -> size_t {
    size_t total = 0;
    for (const auto& buffer : mBufferLists) {
        total += buffer->getPressureEvents();
    }
    return total;
}

void ReadPool::printPressureStats() const {
    size_t totalPressure = getTotalPressureEvents();
    size_t totalProduced = mTotalProduced.load();
    
    if (totalProduced > 0) {
        double pressureRatio = (double)totalPressure / totalProduced * 100.0;
        printf("ReadPool Adaptive Backoff Stats:\n");
        printf("  Pressure events: %zu / %zu (%.2f%%)\n", 
               totalPressure, totalProduced, pressureRatio);
        
        if (pressureRatio > 5.0) {
            printf("High pressure detected. Consider increasing thread count or buffer sizes.\n");
        } else if (pressureRatio < 0.1) {
            printf("Low pressure - adaptive backoff working well!\n");
        } else {
            printf("Moderate pressure - system working efficiently.\n");
        }
    }
}

void ReadPool::resetPressureStats() {
    for (auto& buffer : mBufferLists) {
        buffer->resetPressureEvents();
    }
    mTotalProduced.store(0);
    mTotalConsumed.store(0);
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