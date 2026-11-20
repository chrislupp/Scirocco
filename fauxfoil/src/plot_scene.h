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

/// @file plot_scene.h
/// A plot as a list of strokes and labels on a fixed page, independent of
/// how it is drawn (a Qt window, a PostScript hardcopy).

#include <string>
#include <utility>
#include <vector>

namespace fauxfoil {

/// Plot colours. Fore is the foreground: white on screen, black on paper.
enum class Ink {
    Fore,
    Red,
    Blue,
    DarkRed,
    DarkBlue,
    Green,
    Orange,
    Magenta,
    Cyan,
    Yellow
};

/// Line patterns.
enum class Dash { Solid, Dashed, Dotted };

/// A point on the page.
struct Point {
    double x = 0.0;
    double y = 0.0;
};

/// An axis-aligned rectangle on the page.
struct Rect {
    double x0 = 0.0;
    double y0 = 0.0;
    double x1 = 0.0;
    double y1 = 0.0;
    [[nodiscard]] bool Empty() const { return x1 <= x0 || y1 <= y0; }
    [[nodiscard]] bool Contains(const Point& p) const {
        return p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1;
    }
};

/// A polyline, or a filled polygon, optionally clipped to a rectangle.
struct Stroke {
    Ink ink = Ink::Fore;
    Dash dash = Dash::Solid;
    bool fill = false;
    std::vector<Point> points;
    Rect clip;  ///< Drawn only inside this rectangle unless it is empty
};

/// Linear map between page coordinates and a plot's own (user)
/// coordinates: page = origin + scale * user.
struct UserAxes {
    bool set = false;
    double ox = 0.0;
    double sx = 1.0;
    double oy = 0.0;
    double sy = 1.0;
    [[nodiscard]] Point ToUser(const Point& p) const {
        return {(p.x - ox) / sx, (p.y - oy) / sy};
    }
    [[nodiscard]] Point ToPage(const Point& u) const {
        return {ox + sx * u.x, oy + sy * u.y};
    }
};

/// One piece of a label at one script level.
struct TextRun {
    std::string text;  ///< UTF-8
    int level = 0;     ///< 0 normal, -1 subscript, +1 superscript
};

/// A character placed on the page by LayoutLabel.
struct Glyph {
    std::string text;  ///< One UTF-8 character
    double x = 0.0;    ///< Left edge of its cell
    double y = 0.0;    ///< Baseline
    double height = 0.0;
    double width = 1.0;  ///< Horizontal stretch of the character
    bool slant = false;  ///< Italic
};

/// A text label. The text may contain "_{...}" subscripts and "^{...}"
/// superscripts; everything else, including Greek letters, is UTF-8.
struct Label {
    Ink ink = Ink::Fore;
    Point at;              ///< Baseline point selected by @c align
    double height = 10.0;  ///< Cap height
    double align = 0.0;    ///< 0 left, 0.5 centre, 1 right
    std::string text;
    double width = 1.0;  ///< Character width factor
    bool slant = false;  ///< Italic
};

/// A complete plot. Page coordinates have their origin at the lower left;
/// the page is 1245 by 962 units, the size of XFOIL's default window.
struct Scene {
    static constexpr double kWidth = 1245.0;
    static constexpr double kHeight = 962.0;

    std::vector<Stroke> strokes;
    std::vector<Label> labels;

    /// Part of the page to show (zoom); the whole page when empty.
    Rect view;
    /// The plot's own coordinates (x/c and Cp, x and H, ...) on the page.
    UserAxes user;
    /// Airfoil coordinates on the page, where the plot draws the airfoil.
    UserAxes foil;

    /// Add a polyline.
    void Line(std::vector<Point> points, Ink ink = Ink::Fore,
              Dash dash = Dash::Solid);
    /// Add a straight segment.
    void Segment(double x1, double y1, double x2, double y2,
                 Ink ink = Ink::Fore, Dash dash = Dash::Solid);
    /// Add a filled polygon.
    void Fill(std::vector<Point> points, Ink ink = Ink::Fore);
    /// Add a label.
    void Text(double x, double y, double height, std::string text,
              Ink ink = Ink::Fore, double align = 0.0);
    /// Add a plot symbol of the given type (0 square, 1 circle, 2 triangle,
    /// 3 plus, 4 cross, 5 diamond, 6 yield sign, 7 Y, 8 flipped Y, 9 star,
    /// 10 flipped star, 11 hourglass, 12 bowtie, 13 five-point star)
    /// centred on (x, y), @p size across.
    void Symbol(double x, double y, double size, int type, Ink ink = Ink::Fore);
    /// Scale everything about the page origin (plot-object size).
    void Scale(double factor);
};

/// Page units per "inch": the default page is 11 inches wide.
inline constexpr double kUnitsPerInch = Scene::kWidth / 11.0;

/// Character advance as a fraction of the cap height: the fixed pitch that
/// all renderers use, so a label occupies the same space everywhere.
inline constexpr double kAdvance = 0.85;

/// Split label markup into runs.
std::vector<TextRun> ParseLabel(const std::string& text);

/// Place each character of @p label on the page.
std::vector<Glyph> LayoutLabel(const Label& label);

/// Width of @p label on the page.
double LabelWidth(const Label& label);

/// Replace Latin letters by Greek ones as a Symbol font assigns them
/// (a: alpha, b: beta, ...), for "math characters" annotations.
std::string MathCharacters(const std::string& text);

/// Split UTF-8 text into characters.
std::vector<std::string> Utf8Characters(const std::string& text);

}  // namespace fauxfoil
