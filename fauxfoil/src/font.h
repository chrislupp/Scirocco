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

/// @file font.h
/// fauXFOIL's plot lettering: the FauxfoilPlot TrueType font, built into
/// the program (see fauxfoil/fonts/).

#include <map>
#include <mutex>
#include <string>
#include <string_view>
#include <vector>

#include "plot_scene.h"

namespace fauxfoil {

/// The plot font and the glyph outlines read from it.
class PlotFont {
   public:
    /// The font built into the program.
    static const PlotFont& Get();

    /// Read a TrueType font from @p data, which must outlive the object.
    explicit PlotFont(std::string_view data);

    /// The font file, for registering with a renderer.
    [[nodiscard]] std::string_view Data() const { return data_; }

    /// Cap height as a fraction of the em.
    [[nodiscard]] double CapHeight() const { return cap_height_; }

    /// Whether the font has a glyph for @p code.
    [[nodiscard]] bool Has(char32_t code) const;

    /// Filled outline of a character as closed contours in em units, with
    /// the origin at the left of its cell on the baseline. Characters the
    /// font lacks are drawn as '?'.
    [[nodiscard]] const std::vector<std::vector<Point>>& Outline(
        char32_t code) const;

   private:
    [[nodiscard]] unsigned GlyphIndex(char32_t code) const;
    [[nodiscard]] std::vector<std::vector<Point>> ReadGlyph(
        unsigned index) const;

    std::string_view data_;
    std::map<std::string, std::string_view> tables_;
    double units_per_em_ = 1000.0;
    double cap_height_ = 0.7;
    bool long_loca_ = false;
    unsigned glyph_count_ = 0;
    mutable std::mutex mutex_;
    mutable std::map<char32_t, std::vector<std::vector<Point>>> outlines_;
};

/// The code point of one UTF-8 character.
char32_t CodePoint(const std::string& utf8);

/// Filled outlines of a label's characters on the page.
std::vector<std::vector<Point>> LabelOutlines(const Label& label);

}  // namespace fauxfoil
