#ifndef READ_H_STUB
#define READ_H_STUB
#include <string>

class Read {
public:
    Read(const std::string& name, const std::string& seq, const std::string& qual="")
        : mName(name), mSeq(seq), mQual(qual) {}
    const std::string& name() const { return mName; }
    const std::string& seq() const { return mSeq; }
    int length() const { return static_cast<int>(mSeq.length()); }
private:
    std::string mName;
    std::string mSeq;
    std::string mQual;
};

#endif