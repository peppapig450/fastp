#pragma once

#include <cstddef>
#include <cstring>
#include <string>

namespace fastp {
namespace base {

struct StringRef {
    const char*           data_ = "";
    std::size_t           size_ = 0U;
    static constexpr auto npos  = static_cast<std::size_t>(-1U);

    constexpr StringRef() noexcept = default;

    constexpr StringRef(const char* data, std::size_t size) noexcept
        : data_((data != nullptr) ? data : "")
        , size_(size) {}

    explicit StringRef(const char* c_str) noexcept
        : data_((c_str != nullptr) ? c_str : "")
        , size_((c_str != nullptr) ? std::strlen(c_str) : 0) {}

    // allow from lvalue std::string, forbid rvalue to avoid dangling
    explicit StringRef(const std::string& source) noexcept
        : data_(source.data())
        , size_(source.size()) {}
    explicit StringRef(std::string&&) = delete;

    constexpr auto data() const noexcept -> const char* { return data_; }
    constexpr auto size() const noexcept -> std::size_t { return size_; }
    constexpr auto empty() const noexcept -> bool { return size_ == 0; }

    auto to_string() const -> std::string { return {data_, size_}; }
};
}  // namespace base
}  // namespace fastp
