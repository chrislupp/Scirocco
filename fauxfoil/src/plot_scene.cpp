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
#include "plot_scene.h"

#include <cmath>
#include <map>

namespace fauxfoil {

namespace {

/// Size and baseline shift of sub- and superscripts, relative to the height
/// of the text around them.
constexpr double kScriptScale = 0.7;
constexpr double kSubShift = -0.45;
constexpr double kSuperShift = 0.55;

}  // namespace

void Scene::Line(std::vector<Point> points, Ink ink, Dash dash) {
    if (points.size() < 2) return;
    strokes.push_back({ink, dash, false, std::move(points), Rect{}});
}

void Scene::Segment(double x1, double y1, double x2, double y2, Ink ink,
                    Dash dash) {
    strokes.push_back({ink, dash, false, {{x1, y1}, {x2, y2}}, Rect{}});
}

void Scene::Fill(std::vector<Point> points, Ink ink) {
    if (points.size() < 3) return;
    strokes.push_back({ink, Dash::Solid, true, std::move(points), Rect{}});
}

void Scene::Text(double x, double y, double height, std::string text, Ink ink,
                 double align) {
    labels.push_back({ink, {x, y}, height, align, std::move(text)});
}

void Scene::Symbol(double x, double y, double size, int type, Ink ink) {
    const double r = 0.5 * size;
    auto line = [&](std::initializer_list<Point> unit) {
        std::vector<Point> points;
        for (const Point& p : unit)
            points.push_back({x + r * p.x, y + r * p.y});
        Line(points, ink);
    };
    constexpr double kPi = 3.14159265358979323846;
    switch (type) {
        case 0:
            line({{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {-1, -1}});
            break;
        case 1: {
            std::vector<Point> circle;
            for (int k = 0; k <= 16; ++k) {
                const double t = 2.0 * kPi * k / 16.0;
                circle.push_back({x + r * std::cos(t), y + r * std::sin(t)});
            }
            Line(circle, ink);
            break;
        }
        case 2:
            line({{-1, -0.7}, {1, -0.7}, {0, 1}, {-1, -0.7}});
            break;
        case 3:
            line({{-1, 0}, {1, 0}});
            line({{0, -1}, {0, 1}});
            break;
        case 4:
            line({{-1, -1}, {1, 1}});
            line({{-1, 1}, {1, -1}});
            break;
        case 5:
            line({{0, -1}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}});
            break;
        case 6:
            line({{-1, 0.7}, {1, 0.7}, {0, -1}, {-1, 0.7}});
            break;
        case 7:
            line({{0, -1}, {0, 0}});
            line({{-0.9, 0.9}, {0, 0}, {0.9, 0.9}});
            break;
        case 8:
            line({{0, 1}, {0, 0}});
            line({{-0.9, -0.9}, {0, 0}, {0.9, -0.9}});
            break;
        case 9:
        case 10: {
            const double turn = type == 9 ? 0.0 : kPi / 6.0;
            for (int k = 0; k < 3; ++k) {
                const double t = turn + kPi * k / 3.0;
                line(
                    {{std::cos(t), std::sin(t)}, {-std::cos(t), -std::sin(t)}});
            }
            break;
        }
        case 11:
            line({{-1, 1}, {1, 1}, {-1, -1}, {1, -1}, {-1, 1}});
            break;
        case 12:
            line({{-1, -1}, {-1, 1}, {1, -1}, {1, 1}, {-1, -1}});
            break;
        default: {
            std::vector<Point> star;
            for (int k = 0; k <= 10; ++k) {
                const double t = kPi / 2.0 + kPi * k / 5.0;
                const double rr = k % 2 == 0 ? r : 0.4 * r;
                star.push_back({x + rr * std::cos(t), y + rr * std::sin(t)});
            }
            Line(star, ink);
            break;
        }
    }
}

void Scene::Scale(double factor) {
    if (factor == 1.0) return;
    for (Stroke& stroke : strokes) {
        for (Point& p : stroke.points) {
            p.x *= factor;
            p.y *= factor;
        }
        stroke.clip = {stroke.clip.x0 * factor, stroke.clip.y0 * factor,
                       stroke.clip.x1 * factor, stroke.clip.y1 * factor};
    }
    for (Label& label : labels) {
        label.at.x *= factor;
        label.at.y *= factor;
        label.height *= factor;
    }
    for (UserAxes* axes : {&user, &foil}) {
        axes->ox *= factor;
        axes->oy *= factor;
        axes->sx *= factor;
        axes->sy *= factor;
    }
}

std::vector<std::string> Utf8Characters(const std::string& text) {
    std::vector<std::string> characters;
    for (std::size_t i = 0; i < text.size();) {
        const auto lead = static_cast<unsigned char>(text[i]);
        std::size_t length = 1;
        if (lead >= 0xF0)
            length = 4;
        else if (lead >= 0xE0)
            length = 3;
        else if (lead >= 0xC0)
            length = 2;
        characters.push_back(text.substr(i, length));
        i += length;
    }
    return characters;
}

std::vector<TextRun> ParseLabel(const std::string& text) {
    std::vector<TextRun> runs;
    TextRun current;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if ((c == '_' || c == '^') && i + 1 < text.size() &&
            text[i + 1] == '{') {
            const std::size_t close = text.find('}', i + 2);
            if (close != std::string::npos) {
                if (!current.text.empty()) runs.push_back(current);
                runs.push_back(
                    {text.substr(i + 2, close - i - 2), c == '_' ? -1 : 1});
                current = TextRun{};
                i = close;
                continue;
            }
        }
        current.text += c;
    }
    if (!current.text.empty()) runs.push_back(current);
    return runs;
}

std::vector<Glyph> LayoutLabel(const Label& label) {
    std::vector<Glyph> glyphs;
    double x = 0.0;
    for (const TextRun& run : ParseLabel(label.text)) {
        const double height =
            run.level == 0 ? label.height : kScriptScale * label.height;
        const double shift = run.level < 0   ? kSubShift * label.height
                             : run.level > 0 ? kSuperShift * label.height
                                             : 0.0;
        for (const std::string& character : Utf8Characters(run.text)) {
            glyphs.push_back(
                {character, x, shift, height, label.width, label.slant});
            x += kAdvance * height * label.width;
        }
    }
    const double left = label.at.x - label.align * x;
    for (Glyph& glyph : glyphs) {
        glyph.x += left;
        glyph.y += label.at.y;
    }
    return glyphs;
}

double LabelWidth(const Label& label) {
    double width = 0.0;
    for (const TextRun& run : ParseLabel(label.text))
        width += kAdvance * (run.level == 0 ? 1.0 : kScriptScale) *
                 label.height * label.width *
                 static_cast<double>(Utf8Characters(run.text).size());
    return width;
}

std::string MathCharacters(const std::string& text) {
    /* The usual Symbol-font assignment of Greek letters to Latin keys */
    static const std::map<char, const char*> kGreek = {
        {'a', "α"}, {'b', "β"}, {'c', "χ"}, {'d', "δ"}, {'e', "ε"}, {'f', "φ"},
        {'g', "γ"}, {'h', "η"}, {'i', "ι"}, {'j', "φ"}, {'k', "κ"}, {'l', "λ"},
        {'m', "μ"}, {'n', "ν"}, {'o', "o"}, {'p', "π"}, {'q', "θ"}, {'r', "ρ"},
        {'s', "σ"}, {'t', "τ"}, {'u', "υ"}, {'v', "ν"}, {'w', "ω"}, {'x', "ξ"},
        {'y', "ψ"}, {'z', "ζ"}, {'A', "A"}, {'B', "B"}, {'C', "X"}, {'D', "Δ"},
        {'E', "E"}, {'F', "Φ"}, {'G', "Γ"}, {'H', "H"}, {'I', "I"}, {'J', "Θ"},
        {'K', "K"}, {'L', "Λ"}, {'M', "M"}, {'N', "N"}, {'O', "O"}, {'P', "Π"},
        {'Q', "Θ"}, {'R', "P"}, {'S', "Σ"}, {'T', "T"}, {'U', "Y"}, {'V', "V"},
        {'W', "Ω"}, {'X', "Ξ"}, {'Y', "Ψ"}, {'Z', "Z"},
    };
    std::string out;
    for (char c : text) {
        const auto found = kGreek.find(c);
        out += found != kGreek.end() ? found->second : std::string(1, c);
    }
    return out;
}

}  // namespace fauxfoil
