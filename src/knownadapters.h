#pragma once

#include <string>
#include <unordered_map>

namespace adapters {  // adapters

auto getKnown() -> const std::unordered_map<std::string, std::string>&;
}  // namespace adapters
