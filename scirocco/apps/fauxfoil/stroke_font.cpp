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
#include "stroke_font.h"

#include <cmath>
#include <cstdlib>
#include <map>
#include <sstream>

namespace fauxfoil {

namespace {

/* Glyph outlines. Strokes are separated by ';' and points by blanks. A
   point marked 'r' (or 'R') becomes a rounded corner of radius 1.1 (or
   1.8) glyph units. */
const std::map<std::string, const char*>& Outlines() {
    static const std::map<std::string, const char*> kOutlines = {
        /* Digits */
        {"0", "0,3 0,6r 4,6r 4,0r 0,0r 0,3"},
        {"1", "0.8,4.8 2,6 2,0"},
        {"2", "0,5 0,6r 4,6r 4,4r 0,1.4 0,0 4,0"},
        {"3", "0,5 0,6r 4,6r 4,3.1r 1.4,3.1; 1.4,3.1 4,3.1r 4,0r 0,0r 0,1"},
        {"4", "0,6 0,2.2 4,2.2; 3,6 3,0"},
        {"5", "4,6 0,6 0,3.5 4,3.5r 4,0r 0,0r 0,1"},
        {"6", "4,5 4,6r 0,6r 0,0r 4,0r 4,3.5r 0,3.5"},
        {"7", "0,6 4,6 1.4,0"},
        {"8",
         "2,3.2 0.3,3.2r 0.3,6r 3.7,6r 3.7,3.2r 2,3.2; "
         "2,3.2 0,3.2r 0,0r 4,0r 4,3.2r 2,3.2"},
        {"9", "0,1 0,0r 4,0r 4,6r 0,6r 0,2.5r 4,2.5"},

        /* Capitals */
        {"A", "0,0 0,4 2,6 4,4 4,0; 0,2.6 4,2.6"},
        {"B", "0,0 0,6 3.4,6r 3.4,3.2r 0,3.2; 3,3.2 4,3.2r 4,0r 0,0"},
        {"C", "4,6 0,6r 0,0r 4,0"},
        {"D", "0,0 0,6 4,6R 4,0R 0,0"},
        {"E", "4,6 0,6 0,0 4,0; 0,3.2 3,3.2"},
        {"F", "4,6 0,6 0,0; 0,3.2 3,3.2"},
        {"G", "4,5 4,6r 0,6r 0,0r 4,0r 4,3 2,3"},
        {"H", "0,0 0,6; 4,0 4,6; 0,3.2 4,3.2"},
        {"I", "1,6 3,6; 2,6 2,0; 1,0 3,0"},
        {"J", "4,6 4,0r 0,0r 0,1.5"},
        {"K", "0,0 0,6; 4,6 0,2.2; 1.3,3.4 4,0"},
        {"L", "0,6 0,0 4,0"},
        {"M", "0,0 0,6 2,3 4,6 4,0"},
        {"N", "0,0 0,6 4,0 4,6"},
        {"O", "0,3 0,6R 4,6R 4,0R 0,0R 0,3"},
        {"P", "0,0 0,6 4,6r 4,3r 0,3"},
        {"Q", "0,3 0,6R 4,6R 4,0R 0,0R 0,3; 2.4,1.6 4.2,-0.4"},
        {"R", "0,0 0,6 4,6r 4,3.2r 0,3.2; 2,3.2 4,0"},
        {"S", "4,5 4,6r 0,6r 0,3.2r 4,3.2r 4,0r 0,0r 0,1"},
        {"T", "0,6 4,6; 2,6 2,0"},
        {"U", "0,6 0,0r 4,0r 4,6"},
        {"V", "0,6 2,0 4,6"},
        {"W", "0,6 1,0 2,4 3,0 4,6"},
        {"X", "0,6 4,0; 0,0 4,6"},
        {"Y", "0,6 2,3 4,6; 2,3 2,0"},
        {"Z", "0,6 4,6 0,0 4,0"},

        /* Small letters */
        {"a", "0.4,4 4,4r 4,0; 4,2.2 0,2.2r 0,0r 4,0"},
        {"b", "0,6 0,0 4,0r 4,4r 0,4"},
        {"c", "4,4 0,4r 0,0r 4,0"},
        {"d", "4,6 4,0 0,0r 0,4r 4,4"},
        {"e", "0,2.2 4,2.2 4,4r 0,4r 0,0r 4,0"},
        {"f", "4,6 1.6,6r 1.6,0; 0,4 3.6,4"},
        {"g", "4,0 0,0r 0,4r 4,4 4,-2r 0.4,-2"},
        {"h", "0,6 0,0; 0,4 4,4r 4,0"},
        {"i", "2,4 2,0; 2,5.3 2,5.8"},
        {"j", "2.6,4 2.6,-2r 0.6,-2; 2.6,5.3 2.6,5.8"},
        {"k", "0,6 0,0; 4,4 0,1.6; 1.4,2.4 4,0"},
        {"l", "1,6 2,6 2,0; 1,0 3,0"},
        {"m", "0,0 0,4 2,4r 2,0; 2,4 4,4r 4,0"},
        {"n", "0,0 0,4 4,4r 4,0"},
        {"o", "0,2 0,4r 4,4r 4,0r 0,0r 0,2"},
        {"p", "0,-2 0,4 4,4r 4,0r 0,0"},
        {"q", "4,-2 4,4 0,4r 0,0r 4,0"},
        {"r", "0,0 0,4; 0,2.4 1.6,4r 4,4"},
        {"s", "4,4 0,4r 0,2r 4,2r 4,0r 0,0"},
        {"t", "1.6,6 1.6,0r 4,0; 0,4 3.6,4"},
        {"u", "0,4 0,0r 4,0; 4,4 4,0"},
        {"v", "0,4 2,0 4,4"},
        {"w", "0,4 1,0 2,3 3,0 4,4"},
        {"x", "0,4 4,0; 0,0 4,4"},
        {"y", "0,4 2,0; 4,4 1.2,-2"},
        {"z", "0,4 4,4 0,0 4,0"},

        /* Punctuation and signs */
        {".", "1.6,0 2.4,0 2.4,0.8 1.6,0.8 1.6,0"},
        {",", "2.4,0.8 2.4,0 1.4,-1.2"},
        {":",
         "1.6,0 2.4,0 2.4,0.8 1.6,0.8 1.6,0; "
         "1.6,3.4 2.4,3.4 2.4,4.2 1.6,4.2 1.6,3.4"},
        {";",
         "1.6,3.4 2.4,3.4 2.4,4.2 1.6,4.2 1.6,3.4; 2.4,0.8 2.4,0 1.4,-1.2"},
        {"-", "0.5,3 3.5,3"},
        {"+", "0.5,3 3.5,3; 2,1.5 2,4.5"},
        {"=", "0.8,2.3 3.2,2.3; 0.8,3.7 3.2,3.7"},
        {"/", "0,0 4,6"},
        {"\\", "0,6 4,0"},
        {"(", "3,7 1.5,5r 1.5,1r 3,-1"},
        {")", "1,7 2.5,5r 2.5,1r 1,-1"},
        {"[", "3,7 1.4,7 1.4,-1 3,-1"},
        {"]", "1,7 2.6,7 2.6,-1 1,-1"},
        {"*", "2,1.5 2,4.5; 0.7,2.2 3.3,3.8; 0.7,3.8 3.3,2.2"},
        {"?", "0,5 0,6r 4,6r 4,3.5r 2,3.5 2,2; 2,0 2,0.6"},
        {"!", "2,6 2,1.6; 2,0 2,0.6"},
        {"_", "0,-1 4,-1"},
        {"'", "2,6 2,4.6"},
        {"\"", "1.3,6 1.3,4.6; 2.7,6 2.7,4.6"},
        {"<", "4,5 0,3 4,1"},
        {">", "0,5 4,3 0,1"},
        {"|", "2,7 2,-1"},
        {"#", "1.2,0.5 1.8,5.5; 2.6,0.5 3.2,5.5; 0.2,2 4,2; 0,4 3.8,4"},
        {"%",
         "0,0 4,6; 0.4,4.6 1.4,4.6 1.4,5.6 0.4,5.6 0.4,4.6; "
         "2.6,0.4 3.6,0.4 3.6,1.4 2.6,1.4 2.6,0.4"},
        {"&", "4,0 0.6,4.4 0.6,6r 2.8,6r 2.8,4.4 0,2 0,0r 2.6,0r 4,1.6"},
        {"$",
         "4,4.6 4,5.2r 0,5.2r 0,3.2r 4,3.2r 4,0.8r 0,0.8r 0,1.4; "
         "2,6.2 2,-0.2"},
        {"@", "3,2 3,4 1,4r 1,2r 4,2 4,6r 0,6r 0,0r 4,0"},
        {"^", "0.6,4.6 2,6 3.4,4.6"},
        {"~", "0,3 1,3.8 3,2.4 4,3.2"},
        {"`", "1.4,6.2 2.4,5"},
        {"{", "3,7 2,6.6 2,3.6 1.2,3 2,2.4 2,-0.6 3,-1"},
        {"}", "1,7 2,6.6 2,3.6 2.8,3 2,2.4 2,-0.6 1,-1"},

        /* Greek letters and symbols */
        {"°", "1,5.3 1,6.2r 2.4,6.2r 2.4,4.8r 1,4.8r 1,5.3"},
        {"×", "1.2,2.2 2.8,3.8; 1.2,3.8 2.8,2.2"},
        {"α",
         "4.2,4 3.6,2.2 3,0.9 2.2,0.1 1.2,0 0.4,0.5 0,1.7 0.2,3 "
         "0.9,3.8 1.9,4 2.7,3.2 3.2,1.9 3.7,0.5 4.3,0"},
        {"θ",
         "2,0 0.8,0.4 0.2,1.6 0,3 0.2,4.4 0.8,5.6 2,6 3.2,5.6 3.8,4.4 "
         "4,3 3.8,1.6 3.2,0.4 2,0; 0.1,3.1 3.9,3.1"},
        {"δ",
         "3.6,5.6 2.8,6 1.6,6 1,5.4 1.4,4.8 3.4,3.6 4,2.4 3.8,1 3,0.1 "
         "2,0 1,0.1 0.2,1 0,2.2 0.4,3.4 1.6,4.3"},
        {"τ", "0,3.4 0.8,4 4,4; 2.2,4 2.2,0r 3.4,0"},
        {"ω",
         "0.8,4 0,2.4 0.3,0.6 1,0 1.6,0.4 2,1.8 2.4,0.4 3,0 3.7,0.6 "
         "4,2.4 3.2,4"},
        {"Δ", "0,0 2,6 4,0 0,0"},
        {"∞",
         "2,2 1.4,2.8 0.6,2.8 0,2 0.6,1.2 1.4,1.2 2.6,2.8 3.4,2.8 4,2 "
         "3.4,1.2 2.6,1.2 2,2"},
        {"√", "0,2.4 0.8,2.8 1.8,0 3,6.6 4.6,6.6"},
    };
    return kOutlines;
}

struct Vertex {
    Point p;
    double radius = 0.0;
};

/// Replace marked corners by quadratic arcs.
std::vector<Point> Round(const std::vector<Vertex>& v) {
    std::vector<Point> out;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (v[i].radius <= 0.0 || i == 0 || i + 1 == v.size()) {
            out.push_back(v[i].p);
            continue;
        }
        const Point& a = v[i - 1].p;
        const Point& c = v[i].p;
        const Point& b = v[i + 1].p;
        const double la = std::hypot(a.x - c.x, a.y - c.y);
        const double lb = std::hypot(b.x - c.x, b.y - c.y);
        const double r = std::min({v[i].radius, 0.5 * la, 0.5 * lb});
        if (r <= 0.0) {
            out.push_back(c);
            continue;
        }
        const Point p1{c.x + r * (a.x - c.x) / la, c.y + r * (a.y - c.y) / la};
        const Point p2{c.x + r * (b.x - c.x) / lb, c.y + r * (b.y - c.y) / lb};
        constexpr int kSteps = 5;
        for (int k = 0; k <= kSteps; ++k) {
            const double t = static_cast<double>(k) / kSteps;
            const double u = 1.0 - t;
            out.push_back({u * u * p1.x + 2 * u * t * c.x + t * t * p2.x,
                           u * u * p1.y + 2 * u * t * c.y + t * t * p2.y});
        }
    }
    return out;
}

std::vector<std::vector<Point>> Parse(const char* outline) {
    std::vector<std::vector<Point>> strokes;
    std::stringstream all(outline);
    std::string part;
    while (std::getline(all, part, ';')) {
        std::stringstream points(part);
        std::string token;
        std::vector<Vertex> vertices;
        while (points >> token) {
            Vertex v;
            const char last = token.back();
            if (last == 'r' || last == 'R') {
                v.radius = last == 'r' ? 1.1 : 1.8;
                token.pop_back();
            }
            const std::size_t comma = token.find(',');
            v.p.x = std::strtod(token.substr(0, comma).c_str(), nullptr);
            v.p.y = std::strtod(token.substr(comma + 1).c_str(), nullptr);
            vertices.push_back(v);
        }
        if (vertices.size() > 1) strokes.push_back(Round(vertices));
    }
    return strokes;
}

/// Glyph width in glyph units and units per cap height.
constexpr double kGlyphWidth = 4.0;
constexpr double kCapUnits = 6.0;

}  // namespace

const std::vector<std::vector<Point>>& GlyphStrokes(const std::string& utf8) {
    static std::map<std::string, std::vector<std::vector<Point>>> cache = [] {
        std::map<std::string, std::vector<std::vector<Point>>> parsed;
        for (const auto& [key, outline] : Outlines())
            parsed[key] = Parse(outline);
        parsed[" "] = {};
        return parsed;
    }();
    const auto found = cache.find(utf8);
    if (found != cache.end()) return found->second;
    return cache.at("?");
}

std::vector<std::vector<Point>> LabelStrokes(const Label& label) {
    std::vector<std::vector<Point>> strokes;
    for (const Glyph& glyph : LayoutLabel(label)) {
        const double unit = glyph.height / kCapUnits;
        const double cell = kAdvance * glyph.height;
        const double left = glyph.x + 0.5 * (cell - kGlyphWidth * unit);
        for (const auto& stroke : GlyphStrokes(glyph.text)) {
            std::vector<Point> page;
            page.reserve(stroke.size());
            for (const Point& p : stroke)
                page.push_back({left + unit * p.x, glyph.y + unit * p.y});
            strokes.push_back(std::move(page));
        }
    }
    return strokes;
}

}  // namespace fauxfoil
