#include "fastqreader.h"
#include "read.h"
#include <sys/stat.h>

FastqReader::FastqReader(std::string filename, bool hasQuality, bool phred64) : mIn(filename), mBytesRead(0), mTotalSize(0) {
    mIn.seekg(0, std::ios::end);
    mTotalSize = static_cast<size_t>(mIn.tellg());
    mIn.seekg(0);
}

FastqReader::~FastqReader() {}

Read* FastqReader::read() {
    if(!mIn.good()) return nullptr;
    std::string name, seq, plus, qual;
    if(!std::getline(mIn, name)) return nullptr;
    if(!std::getline(mIn, seq)) return nullptr;
    if(!std::getline(mIn, plus)) return nullptr;
    if(!std::getline(mIn, qual)) return nullptr;
    mBytesRead = static_cast<size_t>(mIn.tellg());
    return new Read(name, seq, qual);
}

void FastqReader::getBytes(size_t& bytesRead, size_t& bytesTotal){
    bytesRead = mBytesRead;
    bytesTotal = mTotalSize;
}

bool FastqReader::eof(){ return mIn.eof(); }
bool FastqReader::hasNoLineBreakAtEnd(){ return false; }
void FastqReader::setReadPool(ReadPool* rp){ }
bool FastqReader::isZipFastq(std::string filename){ return false; }
bool FastqReader::isFastq(std::string filename){ return true; }
bool FastqReader::test(){ return true; }