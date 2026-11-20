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

/// @file format.h
/// Fortran edit descriptors, so fauXFOIL columns line up with XFOIL's output.

#include <string>

namespace fauxfoil {

/// Fortran Fw.d: fixed point, right-justified; '*' fill on overflow.
std::string F(double value, int width, int decimals);

/// Fortran Ew.d: "0.dddE+xx" mantissa form, as gfortran writes it (the
/// leading zero is dropped when the field is too narrow, e.g. "-.7031E+00").
std::string E(double value, int width, int digits);

/// Fortran Gw.d: F editing with four trailing blanks when the magnitude
/// suits it, E editing otherwise.
std::string G(double value, int width, int digits);

/// Fortran Iw: right-justified integer.
std::string I(long value, int width);

/// Fortran Aw on output: left-justified, blank-padded or truncated to width.
std::string A(const std::string& text, int width);

}  // namespace fauxfoil
