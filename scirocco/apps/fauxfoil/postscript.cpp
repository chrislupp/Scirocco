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
#include "postscript.h"

#include <algorithm>
#include <cstdio>
#include <sstream>

#include "stroke_font.h"

namespace fauxfoil {

namespace {

/* Letter paper, landscape, half-inch margins */
constexpr double kPaperLong = 792.0;
constexpr double kPaperShort = 612.0;
constexpr double kMargin = 36.0;

const char* RgbOnPaper(Ink ink) {
    switch (ink) {
        case Ink::Fore:
            return "0 0 0";
        case Ink::Red:
            return "1 0 0";
        case Ink::Blue:
            return "0.118 0.549 1";
        case Ink::DarkRed:
            return "0.667 0 0";
        case Ink::DarkBlue:
            return "0.078 0.365 0.667";
        case Ink::Green:
            return "0 0.6 0";
        case Ink::Orange:
            return "1 0.5 0";
        case Ink::Magenta:
            return "0.8 0 0.8";
        case Ink::Cyan:
            return "0 0.6 0.6";
        case Ink::Yellow:
            return "0.7 0.6 0";
    }
    return "0 0 0";
}

const char* DashPattern(Dash dash) {
    switch (dash) {
        case Dash::Solid:
            return "[] 0";
        case Dash::Dashed:
            return "[3 2] 0";
        case Dash::Dotted:
            return "[0 2.2] 0";
    }
    return "[] 0";
}

}  // namespace

std::string PostScript(const std::vector<Scene>& pages) {
    const double scale =
        std::min((kPaperLong - 2.0 * kMargin) / Scene::kWidth,
                 (kPaperShort - 2.0 * kMargin) / Scene::kHeight);
    const double ox = 0.5 * (kPaperLong - scale * Scene::kWidth);
    const double oy = 0.5 * (kPaperShort - scale * Scene::kHeight);

    std::ostringstream ps;
    ps << "%!PS-Adobe-3.0\n"
       << "%%Title: plot.ps\n"
       << "%%Creator: fauXfoil\n"
       << "%%Pages: " << pages.size() << '\n'
       << "%%Orientation: Landscape\n"
       << "%%BoundingBox: 0 0 612 792\n"
       << "%%EndComments\n"
       << "%%BeginProlog\n"
       << "/M { moveto } bind def\n"
       << "/L { lineto } bind def\n"
       << "%%EndProlog\n";

    char buffer[160];
    int number = 0;
    for (const Scene& scene : pages) {
        ++number;
        ps << "%%Page: " << number << ' ' << number << '\n'
           << "gsave\n"
           << kPaperShort << " 0 translate 90 rotate\n"
           << "1 setlinejoin 1 setlinecap 0.5 setlinewidth\n";
        std::snprintf(buffer, sizeof buffer,
                      "%.3f %.3f translate %.6f %.6f scale\n", ox, oy, scale,
                      scale);
        ps << buffer;
        ps << "/lw " << 0.5 / scale << " def lw setlinewidth\n";

        for (const Stroke& stroke : scene.strokes) {
            ps << RgbOnPaper(stroke.ink) << " setrgbcolor ";
            const std::string dash = DashPattern(stroke.dash);
            if (stroke.dash == Dash::Solid) {
                ps << "[] 0 setdash\n";
            } else {
                /* Dash lengths in points, independent of the page scale */
                std::istringstream parts(dash.substr(1, dash.find(']') - 1));
                ps << '[';
                double v = 0.0;
                while (parts >> v) ps << v / scale << ' ';
                ps << "] 0 setdash\n";
            }
            ps << "newpath\n";
            for (std::size_t i = 0; i < stroke.points.size(); ++i) {
                std::snprintf(buffer, sizeof buffer, "%.2f %.2f %s\n",
                              stroke.points[i].x, stroke.points[i].y,
                              i == 0 ? "M" : "L");
                ps << buffer;
            }
            ps << (stroke.fill ? "closepath fill\n" : "stroke\n");
        }

        ps << "[] 0 setdash\n";
        for (const Label& label : scene.labels) {
            ps << RgbOnPaper(label.ink) << " setrgbcolor\n";
            for (const auto& stroke : LabelStrokes(label)) {
                ps << "newpath";
                for (std::size_t i = 0; i < stroke.size(); ++i) {
                    std::snprintf(buffer, sizeof buffer, " %.2f %.2f %s",
                                  stroke[i].x, stroke[i].y, i == 0 ? "M" : "L");
                    ps << buffer;
                }
                ps << " stroke\n";
            }
        }
        ps << "grestore\nshowpage\n";
    }
    ps << "%%Trailer\n%%EOF\n";
    return ps.str();
}

}  // namespace fauxfoil
