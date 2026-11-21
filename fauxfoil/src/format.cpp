/*
    Scirocco: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

    This work has been cleared for public release, distribution unlimited, case
    number: AFRL-2026-0405.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
#include "format.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace fauxfoil {

namespace {

std::string Stars(int width) {
    std::string stars;
    stars.assign(static_cast<std::size_t>(width), '*');
    return stars;
}

std::string PadLeft(const std::string& text, int width) {
    if (static_cast<int>(text.size()) > width) return Stars(width);
    return std::string(width - text.size(), ' ') + text;
}

}  // namespace

std::string F(double value, int width, int decimals) {
    if (!std::isfinite(value))
        return PadLeft(value != value ? "NaN" : "Inf", width);
    char buffer[128];
    std::snprintf(buffer, sizeof(buffer), "%.*f", decimals, value);
    return PadLeft(buffer, width);
}

std::string E(double value, int width, int digits) {
    if (!std::isfinite(value))
        return PadLeft(value != value ? "NaN" : "Inf", width);
    std::string mantissa;
    int exponent = 0;
    if (value == 0.0) {
        mantissa = std::string(digits, '0');
    } else {
        /* %e gives d.ddd...e+XX with `digits` significant digits; Fortran
           writes the same digits as 0.dddd with the exponent one larger. */
        char buffer[128];
        std::snprintf(buffer, sizeof(buffer), "%.*e", digits - 1,
                      std::fabs(value));
        const std::string text(buffer);
        const std::size_t e = text.find('e');
        for (std::size_t k = 0; k < e; ++k)
            if (text[k] != '.') mantissa += text[k];
        exponent =
            static_cast<int>(std::strtol(text.c_str() + e + 1, nullptr, 10)) +
            1;
    }
    char exp_text[16];
    if (std::abs(exponent) <= 99)
        std::snprintf(exp_text, sizeof(exp_text), "E%c%02d",
                      exponent < 0 ? '-' : '+', std::abs(exponent));
    else
        std::snprintf(exp_text, sizeof(exp_text), "%c%03d",
                      exponent < 0 ? '-' : '+', std::abs(exponent));
    const std::string sign = value < 0.0 ? "-" : "";
    std::string text = sign + "0." + mantissa + exp_text;
    if (static_cast<int>(text.size()) > width)
        text = sign + "." + mantissa + exp_text;
    return PadLeft(text, width);
}

std::string G(double value, int width, int digits) {
    if (value == 0.0) return F(0.0, width - 4, digits - 1) + "    ";
    if (!std::isfinite(value)) return E(value, width, digits);
    char buffer[128];
    std::snprintf(buffer, sizeof(buffer), "%.*e", digits - 1, std::fabs(value));
    const char* e = std::strchr(buffer, 'e');
    const int k = static_cast<int>(std::strtol(e + 1, nullptr, 10)) + 1;
    if (k >= 0 && k <= digits) return F(value, width - 4, digits - k) + "    ";
    return E(value, width, digits);
}

std::string I(long value, int width) {
    return PadLeft(std::to_string(value), width);
}

std::string A(const std::string& text, int width) {
    if (static_cast<int>(text.size()) >= width) return text.substr(0, width);
    return text + std::string(width - text.size(), ' ');
}

}  // namespace fauxfoil
