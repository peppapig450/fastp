#include "readpool.h"
#include "spsc_ring_buffer.h"
#include "util.h"
#include <algorithm>
#include <atomic>
#include <cstddef>
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
    if (data == nullptr) { // nothing to recycle
        return true;
    }

    // try non-blocking produce first
    // This is the most common case and avoids size() calculation
    if (mBufferLists[tid]->tryProduce(data)) {
        // This is for debugging ONLY
        mTotalProduced.fetch_add(1, std::memory_order_relaxed);
        return true;
    }

    // Pool is full, caller should delete the read
    return false;
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
        // This buffer size balanced memory usage and recycling efficiency
        static constexpr size_t kDefaultBufferCapacity = 4UL * 1024 * 1024; // 4M reads per buffer
        size_t bufferCapacity = std::max(kDefaultBufferCapacity, mPerBufferLimit * 2);

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
    size_t totalFailed = 0;

    // Calculate how many reads couldn't be recycled
    for (int tid = 0; tid < mOptions->thread; tid++) {
        size_t bufferSize = mBufferLists[tid]->size();
        if (bufferSize >= mBufferLists[tid]->capacity()) {
            totalFailed++; // Approximate failed recycles
        }
    }
    
    if (totalProduced > 0) {
        double pressureRatio = static_cast<double>(totalPressure) / totalProduced * 100.0;
        printf("ReadPool Stats:\n");
        printf("  Total recycled: %zu\n", totalProduced);
        printf("  Pressure events: %zu / %zu (%.2f%%)\n", 
               totalPressure, totalProduced, pressureRatio);
        
        if (pressureRatio > 5.0) {
            printf("  Note: High pressure detected. Some reads may not be recycled.\n");
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
    // Simple round-robin through all buffers
    int startThread = mLastCheckedThread.load(std::memory_order_relaxed);

    for (int i = 0; i < mOptions->thread; i++) {
        int tid = (startThread + i) % mOptions->thread;

        Read* r = nullptr;
        if (mBufferLists[tid]->tryConsume(r)) {
            mLastCheckedThread.store((tid + 1) % mOptions->thread, std::memory_order_relaxed);
            mTotalConsumed.fetch_add(1, std::memory_order_relaxed);
            return r;
        }
    }

    return nullptr;
}