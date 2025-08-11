#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace benchmark_data {

enum class FragmentModel : std::uint8_t { Normal, Uniform };

// Stringification
constexpr auto to_string(FragmentModel model) noexcept -> std::string_view {
    using enum FragmentModel;

    switch (model) {
            // clang-format off
        case Normal:  return "normal";
        case Uniform: return "uniform";
            // clang-format on
    }
    std::unreachable();
}

// Overload for automatic stringification when inserting into a stream
inline auto operator<<(std::ostream& out_stream, FragmentModel model) -> std::ostream& {
    return out_stream << to_string(model);
}

struct FragmentConfig {
    // defaults
    static constexpr int DEFAULT_STDDEV      = 30;
    static constexpr int DEFAULT_MIN_SIZE    = 100;
    static constexpr int DEFAULT_MEAN        = 300;
    static constexpr int DEFAULT_MAX_SIZE    = 400;
    static constexpr int BUFFER_FOR_MEAN     = 100;
    static constexpr int MINIMUM_READ_LENGTH = 100;
    static constexpr int STDDEV_MULTIPLIER   = 3;
    static constexpr int BUFFER_FOR_MAX_SIZE = 50;

    FragmentModel model = FragmentModel::Normal;

    // For normal distribution
    int mean   = DEFAULT_MEAN;
    int stddev = DEFAULT_STDDEV;

    // For uniform distribution
    int min_size = DEFAULT_MIN_SIZE;
    int max_size = DEFAULT_MAX_SIZE;

    void validate(std::size_t readLen) const {
        const auto intReadLength = static_cast<int>(readLen);

        switch (model) {
            case FragmentModel::Normal:
                if (mean < intReadLength) {
                    throw std::invalid_argument("fragment-mean-size must be >= read length");
                }
                break;

            case FragmentModel::Uniform:
                if (min_size > max_size) {
                    throw std::invalid_argument("fragment-min-size must be <= fragment-max-size");
                }
                if (max_size < intReadLength) {
                    throw std::invalid_argument("fragment-max-size must be >= read length");
                }
                break;
        }
    }

    static constexpr auto SafeDefaults(std::size_t readLen) -> FragmentConfig {
        FragmentConfig config;

        config.model = FragmentModel::Normal;

        const auto intReadLength = static_cast<int>(readLen);

        config.mean     = std::max(intReadLength + BUFFER_FOR_MEAN, DEFAULT_MEAN);
        config.stddev   = DEFAULT_STDDEV;
        config.min_size = std::max(intReadLength, MINIMUM_READ_LENGTH);
        config.max_size = std::max(config.mean + (STDDEV_MULTIPLIER * config.stddev),
                                   intReadLength + BUFFER_FOR_MAX_SIZE);

        return config;
    }
};

}  // namespace benchmark_data
