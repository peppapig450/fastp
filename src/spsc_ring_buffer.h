#pragma once

#include <array>
#include <atomic>
#include <cstddef>
#include <memory>
#include <cassert>
#include <cstdlib>
#include <new>
#include <thread>
#include <type_traits>
#include <chrono>

#if (defined(__x86_64__) || defined(__i386__))                       \
    && (defined(__GNUC__) || defined(__clang__) || defined(_MSC_VER) \
        || defined(__INTEL_LLVM_COMPILER) || defined(__INTEL_COMPILER))
    #include <emmintrin.h>
    #ifndef HAVE__MM_PAUSE
        #define HAVE__MM_PAUSE 1
    #endif
#endif

namespace {
// pipeline-hint to the CPU to pause in a tight spin-loop,
// this avoids stalling the pipeline.
inline void cpu_pause() noexcept {
#if defined(HAVE__MM_PAUSE)
    // x86 PAUSE instruction: hints to CPU this is a spin-wait loop
    // Reduces resource usage and improves hyper-threading performance
    _mm_pause();
#elif defined(__aarch64__) || defined(__arm__)
    // yield instruction for ARM
    asm volatile("yield" ::: "memory");
#elif defined(__riscv)
    // RISC-V single-cycle noop
    asm volatile("nop" ::: "memory");
#elif defined(__powerpc__) || defined(__powerpc64__) || defined(__ppc__) || defined(__PPC__)
    // PowerPC "yield" hint (or 27,27,27): releases shared processor resources
    // see: https://stackoverflow.com/a/7588941
    asm volatile("or 27,27,27" ::: "memory");
#else
    // Fall back for edge cases to a compiler fence, which is not a pause, but prevents
    // reordering.
    asm volatile("" ::: "memory");
    #warning "cpu_pause(): architecture-specific pause not available; using compiler barrier"
#endif
}
}  // namespace

template<typename ItemType>
class SPSCRingBuffer {
private:
    // Cache line size for avoiding false sharing
    static constexpr std::size_t CACHE_LINE_SIZE = 64;

    // Backoff strategy parameters
    static constexpr int MAX_SPINS = 2000;  // Max spins before sleeping
    static constexpr int YIELD_THRESHOLD = 128;   // when to switch from pause to yield
    static constexpr std::chrono::microseconds SLEEP_DURATION{10};

    // Ensure capacity is power of 2 for efficient masking
    static auto roundUpToPowerOf2(std::size_t n) -> std::size_t {
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
        std::atomic<std::size_t> value {0};

        static constexpr std::size_t pad_bytes = CACHE_LINE_SIZE - sizeof(std::atomic<std::size_t>);
        static_assert(pad_bytes > 0, "CACHE_LINE_SIZE too small for std::atomic<std::size_t>");

        // Prevent false sharing with padding
        std::array<char, pad_bytes> padding {};
    };

    void backoff(int& spin_count) {
        backoff_impl(spin_count, [](/* nothing */){});
    }

    template <typename OnLongWaitCallback>
    void backoff(int& spin_count, OnLongWaitCallback&& on_long_wait) {
        backoff_impl(spin_count, std::forward<OnLongWaitCallback>(on_long_wait));
    }

    template <typename OnLongWaitCallback>
    void backoff_impl(int& spin_count, OnLongWaitCallback&& on_long_wait) {
        if (spin_count < YIELD_THRESHOLD) {
            // Phase 1: architecture-specific CPU pause instructions
            // Exponential backoff: more pauses as contention continues
            for (int i = 0; i < (1 << (spin_count / 16)); ++i) {
                cpu_pause();
            }
            ++spin_count;
        } else if (spin_count < MAX_SPINS) {
            // Phase 2: give up time slice but stay ready for quick recovery
            std::this_thread::yield();
            ++spin_count;
        } else {
            // Phase 3: long wait detected, execute callback and sleep
            on_long_wait();             // Execute user callback (inlined)
            std::this_thread::sleep_for(SLEEP_DURATION); // microsecond sleep
            spin_count = 0;             // reset cycle
        }
    }

// TODO: figure out a better tuned default, or have size dynamically calculated by
// calling code.
public:
    explicit SPSCRingBuffer(size_t capacity = static_cast<long>(8) * 1024 * 1024)
        : capacity_(roundUpToPowerOf2(capacity))
        , mask_(capacity_ - 1)
        , buffer_(std::unique_ptr<ItemType[]>(new ItemType[capacity_]))
        , producer_finished_{false}
        , consumer_finished_{false}
        , pressure_events_{0} {

        static_assert(std::is_trivially_move_constructible<ItemType>::value,
                    "T must be trivially move-constructable");

        // Ensure minimum capacity
        assert(capacity_ >= 2);
    }

    ~SPSCRingBuffer() noexcept = default;

    // Ensure allocation is properly aligned for CacheAlignedAtomic
    auto operator new(std::size_t size) -> void* {
        void* ptr = nullptr;
        constexpr std::size_t alignment = alignof(CacheAlignedAtomic);

    #if defined(_MSC_VER)
        ptr = _aligned_malloc(size, alignment);
        if (ptr == nullptr) {
            throw std::bad_alloc();
        }
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

    // Producer interface
    template <typename InputItem>
    void produce(InputItem&& item) {
        int spin_count = 0;
 
        // Rather than always succeeding as the original does which allocates more memory
        // we spin until we can produce (when the buffer has space)
        while (true) {
            const std::size_t head = head_.value.load(std::memory_order_relaxed);
            const std::size_t next_head = head + 1;

            // Check if buffer is full (accounts for keeping 1 slot empty)
            if (next_head - tail_.value.load(std::memory_order_acquire) < capacity_) {
                // Store item using perfect forwarding (zero-copy when possible)
                buffer_[head & mask_] = std::forward<InputItem>(item);

                // Release the item to consumer (ensures item write happens before this)
                head_.value.store(next_head, std::memory_order_release);
                return;
            }

            // Buffer is full: backoff with pressure tracking
            // Lambda captures 'this' to increment pressure counter on long waits
            backoff(spin_count, [this]() {
                pressure_events_.fetch_add(1, std::memory_order_relaxed);
            });
        }
    }

    // Non-blocking version using perfect forwarding
    template <typename InputItem>
    auto tryProduce(InputItem&& item) -> bool {
        const std::size_t head = head_.value.load(std::memory_order_relaxed);
        const std::size_t next_head = head + 1;

        if (next_head - tail_.value.load(std::memory_order_acquire) >= capacity_) {
            return false;
        }

        buffer_[head & mask_] = std::forward<InputItem>(item);
        head_.value.store(next_head, std::memory_order_release);
        return true;
    }

    // Single thread consumer interface
    auto tryConsume(ItemType& item) -> bool {
        const std::size_t tail = tail_.value.load(std::memory_order_relaxed);

        // Check if buffer is empty
        if (tail == head_.value.load(std::memory_order_acquire)) {
            return false; // Buffer empty
        }

        // Load item and release slot to producer
        item = std::move(buffer_[tail & mask_]);
        tail_.value.store(tail +  1, std::memory_order_release);
        return true;
    }

    // Overload of tryConsume for pointer types
    // This version is needed for processor tasks that work directly with pointer types.
    // It avoids ambiguity with the general tryConsume(ItemType&), improves readability,
    // and ensures the buffer contents are correctly moved into the caller's pointer.
    //
    // Enabled only when ItemType is a pointer type.
    template <typename PointerItemType = ItemType>
    auto tryConsume(ItemType*& item) ->
        typename std::enable_if<std::is_pointer<PointerItemType>::value, bool>::type {
        const std::size_t tail = tail_.value.load(std::memory_order_relaxed);

        // Check if buffer is empty
        if (tail == head_.value.load(std::memory_order_acquire)) {
            item = nullptr;
            return false;
        }

        // For pointer types this needs to be handled differently
        // Since we're storing pointers, we can just move the pointer
        item = std::move(buffer_[tail & mask_]);
        tail_.value.store(tail + 1, std::memory_order_release);
        return true;
    }

    // Api compatible blocking consume that returns the item directly
    auto consume() -> ItemType {
        int spin_count = 0;
        while (true) {
            const std::size_t tail = tail_.value.load(std::memory_order_relaxed);

            // Check if buffer has items available
            if (tail != head_.value.load(std::memory_order_acquire)) {
                // Move item out of the buffer
                ItemType item = std::move(buffer_[tail & mask_]);

                // Release the slot to producer
                tail_.value.store(tail + 1, std::memory_order_release);
                return item;
            }

            // Buffer empty: check if producer is done
            if (producer_finished_.load(std::memory_order_acquire)) {
                return ItemType{}; // Return default constructed item if producer finished
            }

            // Buffer emtpy but producer still active: backoff without callback
            // Uses zero-overhead backoff version
            backoff(spin_count);
        }
    }

    // Non-blocking checks
    auto canProduce() const noexcept -> bool {
        const std::size_t head = head_.value.load(std::memory_order_relaxed);
        const std::size_t tail = tail_.value.load(std::memory_order_acquire);
        return (head + 1 - tail) < capacity_;
    }

    auto canConsume() const noexcept -> bool {
        const std::size_t tail = tail_.value.load(std::memory_order_relaxed);
        const std::size_t head = head_.value.load(std::memory_order_acquire);
        return tail != head;
    }

    auto empty() const noexcept -> bool {
        return !canConsume();
    }

    auto size() const noexcept -> std::size_t {
        const std::size_t head = head_.value.load(std::memory_order_acquire);
        const std::size_t tail = tail_.value.load(std::memory_order_acquire);
        return head - tail;
    }

    auto capacity() const noexcept -> std::size_t {
        return capacity_ - 1;   // Keep 1 slot reserved for empty/full distinction
    }

    // Lifecycle management for producer/consumers
    void setProducerFinished() {
        producer_finished_.store(true, std::memory_order_release);
    }

    void setConsumerFinished() {
        consumer_finished_.store(true, std::memory_order_release);
    }

    auto isProducerFinished() const noexcept -> bool {
        return producer_finished_.load(std::memory_order_acquire);
    }

    auto isConsumerFinished() const noexcept -> bool {
        return consumer_finished_.load(std::memory_order_acquire);
    }

    // Alternative interface to match existing API
    auto canBeConsumed() const noexcept -> bool {
        return canConsume();
    }

    // Monitoring and diagnostics
    auto getPressureEvents() const noexcept -> std::size_t {
        return pressure_events_.load(std::memory_order_relaxed);
    }

    void resetPressureEvents() noexcept {
        pressure_events_.store(0, std::memory_order_relaxed);
    }

private:
    const std::size_t capacity_;
    const std::size_t mask_;

    // Buffer storage
    std::unique_ptr<ItemType[]> buffer_;

    // Producer-owned cache line
    alignas(CACHE_LINE_SIZE) CacheAlignedAtomic head_;

    // Consumer-owned cache line
    alignas(CACHE_LINE_SIZE) CacheAlignedAtomic tail_;

    // Status flags
    alignas(CACHE_LINE_SIZE) std::atomic<bool> producer_finished_;
    alignas(CACHE_LINE_SIZE) std::atomic<bool> consumer_finished_;

    // Performance monitoring
    mutable std::atomic_size_t pressure_events_;
};

// Convenience alias
template <typename T>
using SingleProducerSingleConsumerList = SPSCRingBuffer<T>;

// Define the static constexpr member
template <typename ItemType>
constexpr std::chrono::microseconds SPSCRingBuffer<ItemType>::SLEEP_DURATION;