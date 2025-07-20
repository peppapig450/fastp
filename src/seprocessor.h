#ifndef SE_PROCESSOR_H
#define SE_PROCESSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "common.h"
#include "read.h"
#include <cstdlib>
#include <condition_variable>
#include <mutex>
#include <thread>
#include "options.h"
#include "threadconfig.h"
#include "filter.h"
#include "umiprocessor.h"
#include "writerthread.h"
#include "duplicate.h"
#include "spsc_ring_buffer.h"
#include "readpool.h"

using namespace std;

typedef struct ReadRepository ReadRepository;

class SingleEndProcessor{
public:
    SingleEndProcessor(Options* opt);
    ~SingleEndProcessor();
    bool process();

private:
    bool processSingleEnd(ReadPack* pack, ThreadConfig* config);
    void readerTask();
    void processorTask(ThreadConfig* config);
    void initConfig(ThreadConfig* config);
    void initOutput();
    void closeOutput();
    void writerTask(WriterThread* config);
    void recycleToPool(int tid, Read* r);

    static constexpr int InputQueueCapacity = PACK_IN_MEM_LIMIT + 1;

private:
    Options* mOptions;
    atomic_bool mReaderFinished;
    atomic_int mFinishedThreads;
    Filter* mFilter;
    UmiProcessor* mUmiProcessor;
    WriterThread* mLeftWriter;
    WriterThread* mFailedWriter;
    Duplicate* mDuplicate;
    SingleProducerSingleConsumerList<ReadPack*>** mInputLists;
    ReadPool* mReadPool;
};


#endif