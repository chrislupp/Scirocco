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

/// @file stroke_font.h
/// fauXfoil's plot lettering: a single-stroke vector font drawn with the
/// same thin pen as the plot lines, in the style of plotter text.

#include <string>
#include <vector>

#include "plot_scene.h"

namespace fauxfoil {

/// Polylines of one character in glyph units: the glyph occupies x in
/// [0, 4] with the baseline at y = 0, cap height 6, x-height 4 and
/// descenders to y = -2. Unknown characters are drawn as '?'; a space has
/// no strokes.
const std::vector<std::vector<Point>>& GlyphStrokes(const std::string& utf8);

/// Replace Latin letters by Greek ones as a Symbol font assigns them
/// (a: alpha, b: beta, ...), for "math characters" annotations.
std::string MathCharacters(const std::string& text);

/// Polylines of a whole label on the page.
std::vector<std::vector<Point>> LabelStrokes(const Label& label);

}  // namespace fauxfoil
