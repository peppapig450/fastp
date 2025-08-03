#pragma once

#include <string>
#include <unordered_map>

namespace adapters {  // adapters

auto getKnown() -> const std::unordered_map<std::string, std::string>&;
auto matchKnown(const std::string& seq) -> std::string;
}  // namespace adapters
