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
#include "polar.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include "format.h"

namespace fauxfoil {

namespace {

const char* const kColumns =
    "   alpha    CL        CD       CDp       CM     Top_Xtr  Bot_Xtr  "
    "Top_Itr  Bot_Itr";
const char* const kDashes =
    "  ------ -------- --------- --------- -------- -------- -------- "
    "-------- --------";

/// Value after "key =" in a header line, or @p fallback if absent.
double ValueAfter(const std::string& line, const std::string& key,
                  double fallback) {
    const std::size_t at = line.find(key);
    if (at == std::string::npos) return fallback;
    const std::size_t eq = line.find('=', at);
    if (eq == std::string::npos) return fallback;
    const char* begin = line.c_str() + eq + 1;
    char* end = nullptr;
    const double value = std::strtod(begin, &end);
    return end == begin ? fallback : value;
}

std::string Trimmed(const std::string& text) {
    const std::size_t begin = text.find_first_not_of(' ');
    if (begin == std::string::npos) return "";
    const std::size_t end = text.find_last_not_of(' ');
    return text.substr(begin, end - begin + 1);
}

}  // namespace

std::string FormatPolar(const Polar& polar, const std::string& program,
                        const std::string& version) {
    std::ostringstream out;
    out << "  \n";
    out << "       " << A(program, 14) << "Version " << version << '\n';
    out << "  \n";
    out << " Calculated polar for: " << A(polar.name, 48) << '\n';
    out << "  \n";
    out << ' ' << polar.re_type << ' ' << polar.mach_type
        << " Reynolds number fixed          Mach number fixed         \n";
    out << "  \n";
    out << " xtrf = " << F(polar.xtrip_top, 7, 3) << " (top)  "
        << F(polar.xtrip_bot, 11, 3) << " (bottom)  \n";
    out << " Mach = " << F(polar.mach, 7, 3)
        << "     Re = " << F(polar.re / 1.0e6, 9, 3)
        << " e 6     Ncrit =" << F(polar.ncrit_top, 8, 3)
        << F(polar.ncrit_bot, 7, 3) << '\n';
    out << "  \n";
    out << kColumns << '\n';
    out << kDashes << '\n';
    for (const PolarPoint& p : polar.points) {
        out << F(p.alpha, 8, 3) << F(p.cl, 9, 4) << F(p.cd, 10, 5)
            << F(p.cdp, 10, 5) << F(p.cm, 9, 4) << F(p.xtr_top, 9, 4)
            << F(p.xtr_bot, 9, 4) << F(p.itr_top, 9, 4) << F(p.itr_bot, 9, 4)
            << '\n';
    }
    return out.str();
}

bool ReadPolarFile(const std::string& path, Polar& polar) {
    std::ifstream stream(path);
    if (!stream) return false;
    Polar result;
    bool have_header = false;
    bool in_table = false;
    std::string line;
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (in_table) {
            std::istringstream values(line);
            PolarPoint p;
            if (values >> p.alpha >> p.cl >> p.cd >> p.cdp >> p.cm) {
                values >> p.xtr_top >> p.xtr_bot;
                values >> p.itr_top >> p.itr_bot;
                result.points.push_back(p);
            }
            continue;
        }
        const std::size_t name_at = line.find("Calculated polar for:");
        if (name_at != std::string::npos) {
            result.name = Trimmed(line.substr(name_at + 21));
            have_header = true;
        } else if (line.find("xtrf") != std::string::npos) {
            std::istringstream values(line.substr(line.find('=') + 1));
            std::string skip;
            values >> result.xtrip_top >> skip >> result.xtrip_bot;
        } else if (line.find("Mach =") != std::string::npos) {
            result.mach = ValueAfter(line, "Mach", 0.0);
            /* "Re =     1.000 e 6" */
            const std::size_t re_at = line.find("Re =");
            if (re_at != std::string::npos) {
                std::istringstream values(line.substr(re_at + 4));
                double mantissa = 0.0;
                std::string e;
                int exponent = 0;
                values >> mantissa >> e >> exponent;
                result.re = mantissa * std::pow(10.0, exponent);
            }
            const std::size_t nc_at = line.find("Ncrit =");
            if (nc_at != std::string::npos) {
                std::istringstream values(line.substr(nc_at + 7));
                values >> result.ncrit_top;
                if (!(values >> result.ncrit_bot))
                    result.ncrit_bot = result.ncrit_top;
            }
        } else if (line.find("Reynolds number") != std::string::npos) {
            std::istringstream values(line);
            values >> result.re_type >> result.mach_type;
        } else if (line.find("------") != std::string::npos) {
            in_table = true;
        }
    }
    if (!have_header) return false;
    result.save_file = path;
    polar = result;
    return true;
}

}  // namespace fauxfoil
