#pragma once

#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <system_error>
#include <vector>

namespace std {
namespace filesystem {
inline path unique_path(const std::string& model = "%%%%-%%%%-%%%%-%%%%") {
    std::random_device rd;
    std::mt19937       gen(rd());
    const char         hex[] = "0123456789abcdef";
    std::string        s     = model;
    for (char& c : s) {
        if (c == '%') c = hex[gen() & 15];
    }
    return path{s};
}
} // namespace filesystem
} // namespace std

class ScopedFastq {
    std::string _path;

public:
    ScopedFastq(const std::vector<std::string>& names,
                const std::vector<std::string>& seqs) {
        auto tmp = std::filesystem::temp_directory_path() /
                   std::filesystem::unique_path("fastq-%%%%-%%%%-%%%%.fq");
        _path = tmp.string();
        std::ofstream ofs(_path);
        for (size_t i = 0; i < names.size(); ++i) {
            ofs << names[i] << "\n";
            ofs << seqs[i] << "\n";
            ofs << "+\n";
            ofs << std::string(seqs[i].size(), 'I') << "\n";
        }
    }

    const std::string& path() const { return _path; }

    ~ScopedFastq() {
        std::error_code ec;
        std::filesystem::remove(_path, ec);
    }
};
