#pragma once

#include <memory>
#include <string>
#include "writer.h"
#include "options.h"
#include <atomic>
#include "spsc_ring_buffer.h"


class WriterThread{
public:
    WriterThread(Options* opt, const std::string& filename, bool isSTDOUT = false);
    ~WriterThread();

    // disabling copy as we do not two threads managing the same buffer
    WriterThread(const WriterThread&)                    = delete;
    auto operator=(const WriterThread&) -> WriterThread& = delete;

    // disable move as it can lead to half-moved atomics
    WriterThread(WriterThread&&)                    = delete;
    auto operator=(WriterThread&&) -> WriterThread& = delete;

    void initWriter(const std::string& filename1, bool isSTDOUT = false);
    void initBufferLists();

    void cleanup();

    auto isCompleted() const -> bool;
    void output();
    void input(int tid, std::string* data);
    bool setInputCompleted();

    auto getFilename() -> const std::string& { return mFilename; }

private:
    void deleteWriter();
    void advanceWorkingBuffer();

private:
    std::unique_ptr<Writer> mWriter1;
    Options* mOptions;
    std::string mFilename;

    // for spliting output
    std::atomic_bool mInputCompleted{false};
    SingleProducerSingleConsumerList<std::string*>** mBufferLists;
    int mWorkingBufferList;

    static constexpr int BufferListCapacity = PACK_IN_MEM_LIMIT + 1;
};
