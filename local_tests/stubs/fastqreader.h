#ifndef FASTQREADER_STUB_H
#define FASTQREADER_STUB_H
#include <string>
#include <fstream>
#include <cstddef>
class Read;
class ReadPool;
class FastqReader{
public:
    FastqReader(std::string filename, bool hasQuality=true, bool phred64=false);
    ~FastqReader();
    Read* read();
    void getBytes(size_t& bytesRead, size_t& bytesTotal);
    bool eof();
    bool hasNoLineBreakAtEnd();
    void setReadPool(ReadPool* rp);
    static bool isZipFastq(std::string filename);
    static bool isFastq(std::string filename);
    static bool test();
private:
    std::ifstream mIn;
    size_t mBytesRead;
    size_t mTotalSize;
};
#endif