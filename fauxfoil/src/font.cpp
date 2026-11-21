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
#include "font.h"

#include <cstdint>
#include <stdexcept>

#include "font_data.h"

namespace fauxfoil {

namespace {

/// Big-endian reads from a font table.
std::uint32_t U32(std::string_view d, std::size_t at) {
    if (at + 4 > d.size()) throw std::out_of_range("font table");
    return (std::uint32_t{static_cast<unsigned char>(d[at])} << 24U) |
           (std::uint32_t{static_cast<unsigned char>(d[at + 1])} << 16U) |
           (std::uint32_t{static_cast<unsigned char>(d[at + 2])} << 8U) |
           std::uint32_t{static_cast<unsigned char>(d[at + 3])};
}

std::uint16_t U16(std::string_view d, std::size_t at) {
    if (at + 2 > d.size()) throw std::out_of_range("font table");
    return static_cast<std::uint16_t>(
        (static_cast<unsigned char>(d[at]) << 8U) |
        static_cast<unsigned char>(d[at + 1]));
}

std::int16_t S16(std::string_view d, std::size_t at) {
    return static_cast<std::int16_t>(U16(d, at));
}

/// TrueType flag bits of glyph points.
constexpr unsigned kOnCurve = 0x01;
constexpr unsigned kXShort = 0x02;
constexpr unsigned kYShort = 0x04;
constexpr unsigned kRepeat = 0x08;
constexpr unsigned kXSame = 0x10;
constexpr unsigned kYSame = 0x20;

/// A contour from TrueType points: quadratic curves are flattened, with an
/// on-curve point implied between two off-curve ones.
std::vector<Point> Flatten(const std::vector<Point>& points,
                           const std::vector<bool>& on) {
    const std::size_t n = points.size();
    std::vector<Point> out;
    if (n == 0) return out;
    /* Start at an on-curve point, or between the first two off-curve */
    std::size_t first = 0;
    while (first < n && !on[first]) ++first;
    Point start = first < n ? points[first]
                            : Point{0.5 * (points[0].x + points[1 % n].x),
                                    0.5 * (points[0].y + points[1 % n].y)};
    if (first == n) first = 0;
    out.push_back(start);
    Point current = start;
    constexpr int kSteps = 4;
    for (std::size_t k = 1; k <= n; ++k) {
        const std::size_t i = (first + k) % n;
        if (on[i]) {
            current = points[i];
            out.push_back(current);
            continue;
        }
        const Point& control = points[i];
        const std::size_t j = (i + 1) % n;
        const Point end = on[j] ? points[j]
                                : Point{0.5 * (control.x + points[j].x),
                                        0.5 * (control.y + points[j].y)};
        for (int s = 1; s <= kSteps; ++s) {
            const double t = static_cast<double>(s) / kSteps;
            const double u = 1.0 - t;
            out.push_back(
                {u * u * current.x + 2 * u * t * control.x + t * t * end.x,
                 u * u * current.y + 2 * u * t * control.y + t * t * end.y});
        }
        current = end;
        if (on[j]) ++k;
    }
    return out;
}

}  // namespace

const PlotFont& PlotFont::Get() {
    static const PlotFont font(std::string_view(
        reinterpret_cast<const char*>(kPlotFontData), kPlotFontSize));
    return font;
}

PlotFont::PlotFont(std::string_view data) : data_(data) {
    const unsigned count = U16(data_, 4);
    for (unsigned k = 0; k < count; ++k) {
        const std::size_t record = 12 + 16 * static_cast<std::size_t>(k);
        const std::string tag(data_.substr(record, 4));
        const std::uint32_t offset = U32(data_, record + 8);
        const std::uint32_t length = U32(data_, record + 12);
        if (offset + static_cast<std::size_t>(length) > data_.size())
            throw std::runtime_error("plot font: truncated table " + tag);
        tables_[tag] = data_.substr(offset, length);
    }
    for (const char* tag : {"head", "maxp", "cmap", "loca", "glyf"})
        if (tables_.count(tag) == 0)
            throw std::runtime_error(std::string("plot font: no ") + tag +
                                     " table");
    const std::string_view head = tables_["head"];
    units_per_em_ = U16(head, 18);
    long_loca_ = S16(head, 50) != 0;
    glyph_count_ = U16(tables_["maxp"], 4);
    const auto os2 = tables_.find("OS/2");
    if (os2 != tables_.end() && U16(os2->second, 0) >= 2 &&
        os2->second.size() >= 90)
        cap_height_ = S16(os2->second, 88) / units_per_em_;
}

unsigned PlotFont::GlyphIndex(char32_t code) const {
    const std::string_view cmap = tables_.at("cmap");
    const unsigned subtables = U16(cmap, 2);
    for (unsigned k = 0; k < subtables; ++k) {
        const std::size_t record = 4 + 8 * static_cast<std::size_t>(k);
        const unsigned platform = U16(cmap, record);
        const unsigned encoding = U16(cmap, record + 2);
        const std::size_t at = U32(cmap, record + 4);
        const bool unicode = platform == 0 || (platform == 3 && encoding == 1);
        if (!unicode || U16(cmap, at) != 4 || code > 0xFFFF) continue;
        /* Format 4: segments of consecutive codes */
        const unsigned segments = U16(cmap, at + 6) / 2;
        const std::size_t ends = at + 14;
        const std::size_t starts = ends + 2 * segments + 2;
        const std::size_t deltas = starts + 2 * segments;
        const std::size_t ranges = deltas + 2 * segments;
        for (unsigned s = 0; s < segments; ++s) {
            if (code > U16(cmap, ends + 2 * s)) continue;
            const unsigned start = U16(cmap, starts + 2 * s);
            if (code < start) return 0;
            const unsigned delta = U16(cmap, deltas + 2 * s);
            const unsigned range = U16(cmap, ranges + 2 * s);
            if (range == 0) return (code + delta) & 0xFFFFU;
            const std::size_t slot =
                ranges + 2 * s + range + 2 * (code - start);
            const unsigned glyph = U16(cmap, slot);
            return glyph == 0 ? 0 : (glyph + delta) & 0xFFFFU;
        }
    }
    return 0;
}

bool PlotFont::Has(char32_t code) const { return GlyphIndex(code) != 0; }

std::vector<std::vector<Point>> PlotFont::ReadGlyph(unsigned index) const {
    std::vector<std::vector<Point>> contours;
    if (index >= glyph_count_) return contours;
    const std::string_view loca = tables_.at("loca");
    const std::string_view glyf = tables_.at("glyf");
    const std::size_t begin = long_loca_
                                  ? U32(loca, 4 * index)
                                  : 2 * std::size_t{U16(loca, 2 * index)};
    const std::size_t end = long_loca_
                                ? U32(loca, 4 * (index + 1))
                                : 2 * std::size_t{U16(loca, 2 * (index + 1))};
    if (end <= begin) return contours; /* empty glyph, such as a space */
    const std::string_view g = glyf.substr(begin, end - begin);
    const int count = S16(g, 0);
    if (count <= 0) return contours; /* composite glyphs are not used */

    std::vector<unsigned> ends(static_cast<std::size_t>(count));
    for (int c = 0; c < count; ++c) ends[c] = U16(g, 10 + 2 * c);
    const std::size_t points = ends.back() + 1;
    std::size_t at = 10 + 2 * static_cast<std::size_t>(count);
    at += 2 + U16(g, at); /* skip instructions */

    std::vector<unsigned> flags;
    flags.reserve(points);
    while (flags.size() < points) {
        const unsigned flag = static_cast<unsigned char>(g.at(at++));
        flags.push_back(flag);
        if ((flag & kRepeat) != 0U) {
            const unsigned repeat = static_cast<unsigned char>(g.at(at++));
            for (unsigned r = 0; r < repeat; ++r) flags.push_back(flag);
        }
    }
    std::vector<Point> xy(points);
    int value = 0;
    for (std::size_t i = 0; i < points; ++i) {
        if ((flags[i] & kXShort) != 0U) {
            const int d = static_cast<unsigned char>(g.at(at++));
            value += (flags[i] & kXSame) != 0U ? d : -d;
        } else if ((flags[i] & kXSame) == 0U) {
            value += S16(g, at);
            at += 2;
        }
        xy[i].x = value / units_per_em_;
    }
    value = 0;
    for (std::size_t i = 0; i < points; ++i) {
        if ((flags[i] & kYShort) != 0U) {
            const int d = static_cast<unsigned char>(g.at(at++));
            value += (flags[i] & kYSame) != 0U ? d : -d;
        } else if ((flags[i] & kYSame) == 0U) {
            value += S16(g, at);
            at += 2;
        }
        xy[i].y = value / units_per_em_;
    }

    std::size_t first = 0;
    for (unsigned last : ends) {
        std::vector<Point> contour(xy.begin() + static_cast<long>(first),
                                   xy.begin() + static_cast<long>(last) + 1);
        std::vector<bool> on;
        on.reserve(contour.size());
        for (std::size_t i = first; i <= last; ++i)
            on.push_back((flags[i] & kOnCurve) != 0U);
        contours.push_back(Flatten(contour, on));
        first = last + 1;
    }
    return contours;
}

const std::vector<std::vector<Point>>& PlotFont::Outline(char32_t code) const {
    const std::scoped_lock lock(mutex_);
    auto found = outlines_.find(code);
    if (found != outlines_.end()) return found->second;
    unsigned index = GlyphIndex(code);
    if (index == 0 && code != U' ') index = GlyphIndex(U'?');
    return outlines_[code] = ReadGlyph(index);
}

char32_t CodePoint(const std::string& utf8) {
    if (utf8.empty()) return 0;
    const auto lead = static_cast<unsigned char>(utf8[0]);
    auto more = [&](std::size_t k) {
        return k < utf8.size()
                   ? static_cast<char32_t>(static_cast<unsigned char>(utf8[k]) &
                                           0x3FU)
                   : char32_t{0};
    };
    if (lead < 0x80) return lead;
    if (lead >= 0xF0)
        return ((lead & 0x07U) << 18U) | (more(1) << 12U) | (more(2) << 6U) |
               more(3);
    if (lead >= 0xE0)
        return ((lead & 0x0FU) << 12U) | (more(1) << 6U) | more(2);
    return ((lead & 0x1FU) << 6U) | more(1);
}

std::vector<std::vector<Point>> LabelOutlines(const Label& label) {
    const PlotFont& font = PlotFont::Get();
    std::vector<std::vector<Point>> outlines;
    for (const Glyph& glyph : LayoutLabel(label)) {
        const double em = glyph.height / font.CapHeight();
        const double wide = em * glyph.width;
        const double shear = glyph.slant ? 0.25 * em : 0.0;
        for (const auto& contour : font.Outline(CodePoint(glyph.text))) {
            std::vector<Point> page;
            page.reserve(contour.size());
            for (const Point& p : contour)
                page.push_back(
                    {glyph.x + wide * p.x + shear * p.y, glyph.y + em * p.y});
            outlines.push_back(std::move(page));
        }
    }
    return outlines;
}

}  // namespace fauxfoil
