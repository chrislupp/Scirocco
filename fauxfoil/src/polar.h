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
#pragma once

/// @file polar.h
/// Stored polars and the XFOIL 6.99 polar save-file layout.

#include <string>
#include <vector>

namespace fauxfoil {

/// One converged operating point, as written to a polar file.
struct PolarPoint {
    double alpha = 0.0;
    double cl = 0.0;
    double cd = 0.0;
    double cdp = 0.0; /**< Surface pressure drag (near field) */
    double cm = 0.0;
    double xtr_top = 0.0; /**< Transition x/c (0 when inviscid) */
    double xtr_bot = 0.0;
    double itr_top = 0.0; /**< Transition as a fractional XFOIL node index */
    double itr_bot = 0.0;
};

/// A polar: fixed parameters plus accumulated points.
struct Polar {
    std::string name;
    double re = 0.0;
    double mach = 0.0;
    double ncrit_top = 9.0;
    double ncrit_bot = 9.0;
    double xtrip_top = 1.0;
    double xtrip_bot = 1.0;
    int re_type = 1;
    int mach_type = 1;
    std::string save_file; /**< Empty when points are not saved */
    std::vector<PolarPoint> points;
    /** Airfoil archived with the polar (APLO, ASET), Scirocco point order */
    std::vector<double> airfoil_x;
    std::vector<double> airfoil_y;
};

/// Full text of a polar save file; @p program and @p version fill the
/// header line where XFOIL writes its own name.
std::string FormatPolar(const Polar& polar, const std::string& program,
                        const std::string& version);

/// Read a polar save file written by fauXFOIL or XFOIL.
/// @return false if the file cannot be opened or has no polar header.
bool ReadPolarFile(const std::string& path, Polar& polar);

}  // namespace fauxfoil
