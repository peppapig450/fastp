#pragma once

#include <atomic>
#include <memory>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <new>
#include <thread>
#include <type_traits>
#include <utility>


template<typename T>
class SPSCRingBuffer {
private:
    // Cache line size for avoiding false sharing
    static constexpr size_t CACHE_LINE_SIZE = 64;

    // Ensure capacity is power of 2 for efficient masking
    static auto roundUpToPowerOf2(size_t n) -> size_t {
        if (n <= 1) { return 2; }
        n--;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        n |= n >> 32;
        return n + 1;
    }

    // Cache-aligned atomic wrapper
    struct alignas(CACHE_LINE_SIZE) CacheAlignedAtomic {
        std::atomic<size_t> value{0};

        // Prevent false sharing with padding
        char padding[CACHE_LINE_SIZE - sizeof(std::atomic<size_t>)];
    };

public:
    explicit SPSCRingBuffer(size_t capacity = 1024 * 1024)
        : capacity_(roundUpToPowerOf2(capacity))
        , mask_(capacity_ - 1)
        , buffer_(std::unique_ptr<T[]>(new T[capacity_]))
        , producer_finished_{false}
        , consumer_finished_{false} {

        static_assert(std::is_trivially_move_constructible<T>::value,
                    "T must be trivially move-constructable");

        // Ensure minimum capacity
        assert(capacity_ >= 2);
    }

    ~SPSCRingBuffer() = default;

    // Ensure allocation is properly aligned for CacheAlignedAtomic
    auto operator new(std::size_t size) -> void* {
        void* ptr = nullptr;
        constexpr std::size_t alignment = alignof(CacheAlignedAtomic);

    #if defined(_MSC_VER)
        ptr = _aligned_malloc(size, alignment);
        if (!ptr) throw std::bad_alloc;
    #else
        if (posix_memalign(&ptr, alignment, size) != 0) {
            throw std::bad_alloc();
        }
    #endif
        return ptr;
    }

    void operator delete(void* ptr) noexcept {
    #if defined(_MSC_VER)
        _aligned_free(ptr);
    #else
        std::free(ptr);
    #endif
    }

    // Non-copyable, non-movable
    SPSCRingBuffer(const SPSCRingBuffer&) = delete;
    auto operator=(const SPSCRingBuffer&) -> SPSCRingBuffer& = delete;
    SPSCRingBuffer(SPSCRingBuffer&&) = delete;
    auto operator=(SPSCRingBuffer&&) -> SPSCRingBuffer& = delete;

    // Producer interface (single thread only)
    void produce(const T& item) {
        // Rather than always succeeding as the original does which allocates more memory
        // we spin until we can produce (when the buffer has space)
        while (true) {
            const size_t head = head_.value.load(std::memory_order_relaxed);
            const size_t next_head = head + 1;

            // Check if buffer is full
            if (next_head - tail_.value.load(std::memory_order_acquire) < capacity_) {
                // Store the item
                buffer_[head & mask_] = item;

                // Release the item to consumer
                head_.value.store(next_head, std::memory_order_release);
                return;
            }

            // Buffer is full, yield and retry
            std::this_thread::yield();
        }
    }

    void produce(T&& item) {
        while (true) {
            const size_t head = head_.value.load(std::memory_order_relaxed);
            const size_t next_head = head + 1;

            if (next_head - tail_.value.load(std::memory_order_acquire) < capacity_) {
                buffer_[head & mask_] = std::move(item);
                head_.value.store(next_head, std::memory_order_release);
                return;
            }

            std::this_thread::yield();
        }
    }

    // Non-blocking versions
    auto tryProduce(const T& item) -> bool {
        const size_t head = head_.value.load(std::memory_order_relaxed);
        const size_t next_head = head + 1;

        if (next_head - tail_.value.load(std::memory_order_acquire) >= capacity_) {
            return false;
        }

        buffer_[head & mask_] = item;
        head_.value.store(next_head, std::memory_order_release);
        return true;
    }

    // Single thread consumer interface
    auto tryConsume(T& item) -> bool {
        const size_t tail = tail_.value.load(std::memory_order_relaxed);

        // Check if buffer is empty
        if (tail == head_.value.load(std::memory_order_acquire)) {
            return false; // Buffer empty
        }

        // Load item and release slot to producer
        item = std::move(buffer_[tail & mask_]);
        tail_.value.store(tail +  1, std::memory_order_release);
        return true;
    }

    // Api compatible blocking consume that returns the item directly
    auto consume() -> T {
        const size_t tail = tail_.value.load(std::memory_order_relaxed);

        // Check if buffer is empty, this should only be called when canBeConsumed() is true
        if (tail == head_.value.load(std::memory_order_acquire)) {
            assert(false && "consume() called on empty buffer");
            return T{};
        }

        T item = std::move(buffer_[tail & mask_]);

        tail_.value.store(tail + 1, std::memory_order_release);
        return item;
    }

    // Non-blocking checks
    auto canProduce() const -> bool {
        const size_t head = head_.value.load(std::memory_order_relaxed);
        const size_t tail = tail_.value.load(std::memory_order_acquire);
        return (head + 1 - tail) < capacity_;
    }

    auto canConsume() const -> bool {
        const size_t tail = tail_.value.load(std::memory_order_relaxed);
        const size_t head = head_.value.load(std::memory_order_acquire);
        return tail != head;
    }

    auto empty() const -> bool {
        const size_t tail = tail_.value.load(std::memory_order_relaxed);
        const size_t head = head_.value.load(std::memory_order_acquire);
        return tail == head;
    }

    auto size() const -> size_t {
        const size_t head = head_.value.load(std::memory_order_acquire);
        const size_t tail = tail_.value.load(std::memory_order_acquire);
        return head - tail;
    }

    auto capacity() const -> size_t {
        return capacity_ - 1;   // Keep 1 slot reserved for empty/full distinction
    }

    // Lifecycle management for producer/consumers
    void setProducerFinished() {
        producer_finished_.store(true, std::memory_order_release);
    }

    void setConsumerFinished() {
        consumer_finished_.store(true, std::memory_order_release);
    }

    auto isProducerFinished() const -> bool {
        return producer_finished_.load(std::memory_order_acquire);
    }

    auto isConsumerFinished() const -> bool {
        return consumer_finished_.load(std::memory_order_acquire);
    }

    // Alternative interface to match existing API
    auto canBeConsumed() const -> bool {
        return canConsume();
    }

private:
    const size_t capacity_;
    const size_t mask_;

    // Buffer storage
    std::unique_ptr<T[]> buffer_;

    // Producer-owned cache line
    alignas(CACHE_LINE_SIZE) CacheAlignedAtomic head_;

    // Consumer-owned cache line
    alignas(CACHE_LINE_SIZE) CacheAlignedAtomic tail_;

    // Status flags
    alignas(CACHE_LINE_SIZE) std::atomic<bool> producer_finished_;
    std::atomic<bool> consumer_finished_;
};

// Convenience alias
template<typename T>
using SingleProducerSingleConsumerList = SPSCRingBuffer<T>;
