#include "writerthread.h"

#include <memory.h>

#include <atomic>

#include "spsc_ring_buffer.h"
#include "writer.h"

// 0 ~ mOptions->thread-1
WriterThread::WriterThread(Options* opt, const std::string& filename, bool isSTDOUT)
    : mOptions(opt), mFilename(filename), mBufferLists(nullptr), mWorkingBufferList(0) {
    initWriter(filename, isSTDOUT);
    initBufferLists();
}

WriterThread::~WriterThread() { cleanup(); }

auto WriterThread::isCompleted() const -> bool {
    if (!mInputCompleted.load(std::memory_order_acquire)) {
        return false;
    }

    if (mBufferLists == nullptr) {
        return true;
    }

    for (int threadIdx = 0; threadIdx < mOptions->thread; ++threadIdx) {
        auto* bufferList = mBufferLists[threadIdx];
        if (bufferList != nullptr && (!bufferList->empty() || !bufferList->isProducerFinished())) {
            return false;
        }
    }
    return true;
}

auto WriterThread::setInputCompleted() -> bool {
    mInputCompleted.store(true, std::memory_order_release);

    if (mBufferLists != nullptr) {
        for (int threadIdx = 0; threadIdx < mOptions->thread; ++threadIdx) {
            auto* bufferList = mBufferLists[threadIdx];
            if (bufferList != nullptr) {
                bufferList->setProducerFinished();
            }
        }
    }
    return true;
}

void WriterThread::output() {
    if (mBufferLists == nullptr || mWriter1 == nullptr) {
        return;
    }

    auto* list = mBufferLists[mWorkingBufferList];
    if (list == nullptr) {
        // This should never happen with proper initialization but we handle it defensively
        // rather than crashing.
        advanceWorkingBuffer();
        return;
    }

    // NOTE: if we populate the buffer with std::unique_ptr's we wouldn't have do this.
    std::string* str = nullptr;
    if (list->tryConsume(str)) {
        mWriter1->write(str->data(), str->length());
        delete str;
    }
    advanceWorkingBuffer();
}

void WriterThread::input(int tid, std::string* data) {
    if (mBufferLists != nullptr && tid >= 0 && tid < mOptions->thread && mBufferLists[tid] != nullptr) {
        mBufferLists[tid]->produce(data);
    }
}

void WriterThread::cleanup() {
    mWriter1.reset();

    if (mBufferLists != nullptr) {
        for (int threadIdx = 0; threadIdx < mOptions->thread; ++threadIdx) {
            delete mBufferLists[threadIdx];
        }
        delete[] mBufferLists;
        mBufferLists = nullptr;
    }
}

void WriterThread::deleteWriter() { mWriter1.reset(); }

void WriterThread::initWriter(const std::string& filename1, bool isSTDOUT) {
    mWriter1.reset(new Writer(mOptions, filename1, mOptions->compression, isSTDOUT));
}

void WriterThread::initBufferLists() {
    if (mBufferLists != nullptr) {
        cleanup();
    }

    auto threadNum = mOptions->thread;
    mBufferLists   = new SingleProducerSingleConsumerList<std::string*>*[threadNum];
    for (int threadIdx = 0; threadIdx < threadNum; ++threadIdx) {
        mBufferLists[threadIdx] = new SingleProducerSingleConsumerList<std::string*>(BufferListCapacity);
    }
}

void WriterThread::advanceWorkingBuffer() { mWorkingBufferList = (mWorkingBufferList + 1) % mOptions->thread; }
