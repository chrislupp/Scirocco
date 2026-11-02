#pragma once

#include <google/protobuf/struct.pb.h>

#include <cmath>
#include <stdexcept>
#include <string>

namespace scirocco {
inline double NumberOption(const google::protobuf::Value& value,
                           const std::string& name) {
    if (!value.has_number_value() || !std::isfinite(value.number_value()))
        throw std::invalid_argument(name + " must be a finite number");
    return value.number_value();
}

inline int IntegerOption(const google::protobuf::Value& value,
                         const std::string& name, int minimum, int maximum) {
    const double number = NumberOption(value, name);
    if (number < minimum || number > maximum || std::floor(number) != number)
        throw std::invalid_argument(name + " must be an integer between " +
                                    std::to_string(minimum) + " and " +
                                    std::to_string(maximum));
    return static_cast<int>(number);
}
}  // namespace scirocco
