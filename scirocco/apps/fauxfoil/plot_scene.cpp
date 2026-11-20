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
    strokes.push_back({ink, dash, false, std::move(points)});
}

void Scene::Segment(double x1, double y1, double x2, double y2, Ink ink,
                    Dash dash) {
    strokes.push_back({ink, dash, false, {{x1, y1}, {x2, y2}}});
}

void Scene::Fill(std::vector<Point> points, Ink ink) {
    if (points.size() < 3) return;
    strokes.push_back({ink, Dash::Solid, true, std::move(points)});
}

void Scene::Text(double x, double y, double height, std::string text, Ink ink,
                 double align) {
    labels.push_back({ink, {x, y}, height, align, std::move(text)});
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
            glyphs.push_back({character, x, shift, height});
            x += kAdvance * height;
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
                 label.height *
                 static_cast<double>(Utf8Characters(run.text).size());
    return width;
}

}  // namespace fauxfoil
