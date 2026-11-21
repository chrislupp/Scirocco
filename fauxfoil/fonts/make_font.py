#!/usr/bin/env python3
# Scirocco: C++ Airfoil Analysis Library
# Copyright 2024-2026 Christopher A. Lupp
# Licensed under the Apache License, Version 2.0; see LICENSE at the
# repository root.
"""Generate FauxfoilPlot.ttf, fauXFOIL's plot lettering.

The letters are single strokes on a small grid, in the style of plotter
text. TrueType holds filled outlines only, so each stroke is drawn with a
round pen of fixed width and the result is stored as the glyph outline.

Each glyph below lists strokes separated by ';', each a run of "x,y"
points: x from 0 to 4, baseline at y = 0, cap height 6, x-height 4,
descenders to -2. A point marked 'r' (or 'R') becomes a rounded corner of
radius 1.1 (or 1.8) grid units.

Regenerate the font after changing a glyph:

    uv run --no-project --with fonttools --with shapely \\
        python fauxfoil/fonts/make_font.py
"""

import math
import pathlib

from fontTools.fontBuilder import FontBuilder
from fontTools.pens.ttGlyphPen import TTGlyphPen
from shapely.geometry import LineString
from shapely.geometry.polygon import orient
from shapely.ops import unary_union

UNITS_PER_EM = 1000
CAP_HEIGHT = 700
GRID = CAP_HEIGHT / 6.0  # font units per grid unit
ADVANCE = round(0.85 * CAP_HEIGHT)  # fixed pitch, as the plots lay text out
LEFT = (ADVANCE - 4.0 * GRID) / 2.0  # centres the 4-unit glyph in its cell
PEN = 44  # stroke width in font units: about one pixel at 16 px cap height
ASCENT = 900
DESCENT = 300

GLYPHS = {
    "0": "0,3 0,6r 4,6r 4,0r 0,0r 0,3",
    "1": "0.8,4.8 2,6 2,0",
    "2": "0,5 0,6r 4,6r 4,4r 0,1.4 0,0 4,0",
    "3": "0,5 0,6r 4,6r 4,3.1r 1.4,3.1; 1.4,3.1 4,3.1r 4,0r 0,0r 0,1",
    "4": "0,6 0,2.2 4,2.2; 3,6 3,0",
    "5": "4,6 0,6 0,3.5 4,3.5r 4,0r 0,0r 0,1",
    "6": "4,5 4,6r 0,6r 0,0r 4,0r 4,3.5r 0,3.5",
    "7": "0,6 4,6 1.4,0",
    "8": "2,3.2 0.3,3.2r 0.3,6r 3.7,6r 3.7,3.2r 2,3.2; 2,3.2 0,3.2r 0,0r 4,0r 4,3.2r 2,3.2",
    "9": "0,1 0,0r 4,0r 4,6r 0,6r 0,2.5r 4,2.5",
    "A": "0,0 0,4 2,6 4,4 4,0; 0,2.6 4,2.6",
    "B": "0,0 0,6 3.4,6r 3.4,3.2r 0,3.2; 3,3.2 4,3.2r 4,0r 0,0",
    "C": "4,6 0,6r 0,0r 4,0",
    "D": "0,0 0,6 4,6R 4,0R 0,0",
    "E": "4,6 0,6 0,0 4,0; 0,3.2 3,3.2",
    "F": "4,6 0,6 0,0; 0,3.2 3,3.2",
    "G": "4,5 4,6r 0,6r 0,0r 4,0r 4,3 2,3",
    "H": "0,0 0,6; 4,0 4,6; 0,3.2 4,3.2",
    "I": "1,6 3,6; 2,6 2,0; 1,0 3,0",
    "J": "4,6 4,0r 0,0r 0,1.5",
    "K": "0,0 0,6; 4,6 0,2.2; 1.3,3.4 4,0",
    "L": "0,6 0,0 4,0",
    "M": "0,0 0,6 2,3 4,6 4,0",
    "N": "0,0 0,6 4,0 4,6",
    "O": "0,3 0,6R 4,6R 4,0R 0,0R 0,3",
    "P": "0,0 0,6 4,6r 4,3r 0,3",
    "Q": "0,3 0,6R 4,6R 4,0R 0,0R 0,3; 2.4,1.6 4.2,-0.4",
    "R": "0,0 0,6 4,6r 4,3.2r 0,3.2; 2,3.2 4,0",
    "S": "4,5 4,6r 0,6r 0,3.2r 4,3.2r 4,0r 0,0r 0,1",
    "T": "0,6 4,6; 2,6 2,0",
    "U": "0,6 0,0r 4,0r 4,6",
    "V": "0,6 2,0 4,6",
    "W": "0,6 1,0 2,4 3,0 4,6",
    "X": "0,6 4,0; 0,0 4,6",
    "Y": "0,6 2,3 4,6; 2,3 2,0",
    "Z": "0,6 4,6 0,0 4,0",
    "a": "0.4,4 4,4r 4,0; 4,2.2 0,2.2r 0,0r 4,0",
    "b": "0,6 0,0 4,0r 4,4r 0,4",
    "c": "4,4 0,4r 0,0r 4,0",
    "d": "4,6 4,0 0,0r 0,4r 4,4",
    "e": "0,2.2 4,2.2 4,4r 0,4r 0,0r 4,0",
    "f": "4,6 1.6,6r 1.6,0; 0,4 3.6,4",
    "g": "4,0 0,0r 0,4r 4,4 4,-2r 0.4,-2",
    "h": "0,6 0,0; 0,4 4,4r 4,0",
    "i": "2,4 2,0; 2,5.3 2,5.8",
    "j": "2.6,4 2.6,-2r 0.6,-2; 2.6,5.3 2.6,5.8",
    "k": "0,6 0,0; 4,4 0,1.6; 1.4,2.4 4,0",
    "l": "1,6 2,6 2,0; 1,0 3,0",
    "m": "0,0 0,4 2,4r 2,0; 2,4 4,4r 4,0",
    "n": "0,0 0,4 4,4r 4,0",
    "o": "0,2 0,4r 4,4r 4,0r 0,0r 0,2",
    "p": "0,-2 0,4 4,4r 4,0r 0,0",
    "q": "4,-2 4,4 0,4r 0,0r 4,0",
    "r": "0,0 0,4; 0,2.4 1.6,4r 4,4",
    "s": "4,4 0,4r 0,2r 4,2r 4,0r 0,0",
    "t": "1.6,6 1.6,0r 4,0; 0,4 3.6,4",
    "u": "0,4 0,0r 4,0; 4,4 4,0",
    "v": "0,4 2,0 4,4",
    "w": "0,4 1,0 2,3 3,0 4,4",
    "x": "0,4 4,0; 0,0 4,4",
    "y": "0,4 2,0; 4,4 1.2,-2",
    "z": "0,4 4,4 0,0 4,0",
    ".": "1.6,0 2.4,0 2.4,0.8 1.6,0.8 1.6,0",
    ",": "2.4,0.8 2.4,0 1.4,-1.2",
    ":": "1.6,0 2.4,0 2.4,0.8 1.6,0.8 1.6,0; 1.6,3.4 2.4,3.4 2.4,4.2 1.6,4.2 1.6,3.4",
    ";": "1.6,3.4 2.4,3.4 2.4,4.2 1.6,4.2 1.6,3.4; 2.4,0.8 2.4,0 1.4,-1.2",
    "-": "0.5,3 3.5,3",
    "+": "0.5,3 3.5,3; 2,1.5 2,4.5",
    "=": "0.8,2.3 3.2,2.3; 0.8,3.7 3.2,3.7",
    "/": "0,0 4,6",
    "\\": "0,6 4,0",
    "(": "3,7 1.5,5r 1.5,1r 3,-1",
    ")": "1,7 2.5,5r 2.5,1r 1,-1",
    "[": "3,7 1.4,7 1.4,-1 3,-1",
    "]": "1,7 2.6,7 2.6,-1 1,-1",
    "*": "2,1.5 2,4.5; 0.7,2.2 3.3,3.8; 0.7,3.8 3.3,2.2",
    "?": "0,5 0,6r 4,6r 4,3.5r 2,3.5 2,2; 2,0 2,0.6",
    "!": "2,6 2,1.6; 2,0 2,0.6",
    "_": "0,-1 4,-1",
    "'": "2,6 2,4.6",
    "\"": "1.3,6 1.3,4.6; 2.7,6 2.7,4.6",
    "<": "4,5 0,3 4,1",
    ">": "0,5 4,3 0,1",
    "|": "2,7 2,-1",
    "#": "1.2,0.5 1.8,5.5; 2.6,0.5 3.2,5.5; 0.2,2 4,2; 0,4 3.8,4",
    "%": "0,0 4,6; 0.4,4.6 1.4,4.6 1.4,5.6 0.4,5.6 0.4,4.6; 2.6,0.4 3.6,0.4 3.6,1.4 2.6,1.4 2.6,0.4",
    "&": "4,0 0.6,4.4 0.6,6r 2.8,6r 2.8,4.4 0,2 0,0r 2.6,0r 4,1.6",
    "$": "4,4.6 4,5.2r 0,5.2r 0,3.2r 4,3.2r 4,0.8r 0,0.8r 0,1.4; 2,6.2 2,-0.2",
    "@": "3,2 3,4 1,4r 1,2r 4,2 4,6r 0,6r 0,0r 4,0",
    "^": "0.6,4.6 2,6 3.4,4.6",
    "~": "0,3 1,3.8 3,2.4 4,3.2",
    "`": "1.4,6.2 2.4,5",
    "{": "3,7 2,6.6 2,3.6 1.2,3 2,2.4 2,-0.6 3,-1",
    "}": "1,7 2,6.6 2,3.6 2.8,3 2,2.4 2,-0.6 1,-1",
    "°": "1,5.3 1,6.2r 2.4,6.2r 2.4,4.8r 1,4.8r 1,5.3",
    "×": "1.2,2.2 2.8,3.8; 1.2,3.8 2.8,2.2",
    "α": "4.2,4 3.6,2.2 3,0.9 2.2,0.1 1.2,0 0.4,0.5 0,1.7 0.2,3 0.9,3.8 1.9,4 2.7,3.2 3.2,1.9 3.7,0.5 4.3,0",
    "θ": "2,0 0.8,0.4 0.2,1.6 0,3 0.2,4.4 0.8,5.6 2,6 3.2,5.6 3.8,4.4 4,3 3.8,1.6 3.2,0.4 2,0; 0.1,3.1 3.9,3.1",
    "δ": "3.6,5.6 2.8,6 1.6,6 1,5.4 1.4,4.8 3.4,3.6 4,2.4 3.8,1 3,0.1 2,0 1,0.1 0.2,1 0,2.2 0.4,3.4 1.6,4.3",
    "τ": "0,3.4 0.8,4 4,4; 2.2,4 2.2,0r 3.4,0",
    "ω": "0.8,4 0,2.4 0.3,0.6 1,0 1.6,0.4 2,1.8 2.4,0.4 3,0 3.7,0.6 4,2.4 3.2,4",
    "Δ": "0,0 2,6 4,0 0,0",
    "∞": "2,2 1.4,2.8 0.6,2.8 0,2 0.6,1.2 1.4,1.2 2.6,2.8 3.4,2.8 4,2 3.4,1.2 2.6,1.2 2,2",
    "√": "0,2.4 0.8,2.8 1.8,0 3,6.6 4.6,6.6",
    "β": "0,-2 0,6 3.4,6r 3.4,3.2r 1,3.2; 1,3.2 4,3.2r 4,0r 0,0",
    "γ": "0,4 2,0.5 4,4; 2,0.5 2,-2",
    "ε": "4,4 0,4r 0,0r 4,0; 0,2 3,2",
    "ζ": "0.5,6 3.5,6 0,1.5 0,0r 3,0r 3,-1.6 2,-2",
    "η": "0,0 0,4 4,4r 4,-2",
    "ι": "1.5,4 1.5,0r 3,0",
    "κ": "0,4 0,0; 4,4 0,1.6; 1.4,2.4 4,0",
    "λ": "0,0 2,4; 0.6,6 4,0",
    "μ": "0,-2 0,4; 0,0 4,0r 4,4",
    "ν": "0,4 2,0 4,3 4,4",
    "ξ": "3.5,6 0.8,6 0.8,4.2 3,4.2; 3,4.2 0,4.2r 0,1.6r 3.6,1.6r 3.6,-2 2.4,-2",
    "π": "0.2,4 4,4; 1.2,4 1.2,0; 3,4 3,0",
    "ρ": "0,-2 0,4r 4,4r 4,0r 0,0",
    "σ": "4,4 0,4r 0,0r 3.2,0r 3.2,4",
    "υ": "0,4 0,0r 4,0r 4,4",
    "φ": "2,6 2,-2; 2,4 0,4r 0,0r 4,0r 4,4r 2,4",
    "χ": "0,4 4,-2; 0,-2 4,4",
    "ψ": "2,5 2,-2; 0,4 0,0r 4,0r 4,4",
    "Γ": "0,0 0,6 4,6",
    "Θ": "0,3 0,6R 4,6R 4,0R 0,0R 0,3; 1,3 3,3",
    "Λ": "0,0 2,6 4,0",
    "Ξ": "0,6 4,6; 0.8,3 3.2,3; 0,0 4,0",
    "Π": "0,0 0,6 4,6 4,0",
    "Σ": "4,6 0,6 2,3 0,0 4,0",
    "Φ": "2,6.5 2,-0.5; 2,5 0,5r 0,1r 4,1r 4,5r 2,5",
    "Ψ": "2,6 2,0; 0,6 0,3r 4,3r 4,6",
    "Ω": "0,0 1.4,0 1.4,0.8 0,2.4 0,6r 4,6r 4,2.4 2.6,0.8 2.6,0 4,0",
}


def rounded(points):
    """Replace marked corners by quadratic arcs, as a plotter pen would."""
    out = []
    for i, (x, y, radius) in enumerate(points):
        if radius <= 0.0 or i == 0 or i == len(points) - 1:
            out.append((x, y))
            continue
        ax, ay, _ = points[i - 1]
        bx, by, _ = points[i + 1]
        la = math.hypot(ax - x, ay - y)
        lb = math.hypot(bx - x, by - y)
        r = min(radius, 0.5 * la, 0.5 * lb)
        if r <= 0.0:
            out.append((x, y))
            continue
        p1 = (x + r * (ax - x) / la, y + r * (ay - y) / la)
        p2 = (x + r * (bx - x) / lb, y + r * (by - y) / lb)
        steps = 5
        for k in range(steps + 1):
            t = k / steps
            u = 1.0 - t
            out.append((u * u * p1[0] + 2 * u * t * x + t * t * p2[0],
                        u * u * p1[1] + 2 * u * t * y + t * t * p2[1]))
    return out


def strokes(outline):
    """Polylines of a glyph, in font units."""
    result = []
    for part in outline.split(";"):
        points = []
        for token in part.split():
            radius = 0.0
            if token[-1] in "rR":
                radius = 1.1 if token[-1] == "r" else 1.8
                token = token[:-1]
            x, y = (float(v) for v in token.split(","))
            points.append((x, y, radius))
        if len(points) > 1:
            result.append([(LEFT + GRID * x, GRID * y)
                           for x, y in rounded(points)])
    return result


def outline_shape(outline):
    """The area a round pen covers along the strokes."""
    pens = [LineString(s).buffer(PEN / 2.0, quad_segs=4, cap_style="round",
                                 join_style="round")
            for s in strokes(outline)]
    return unary_union(pens).simplify(0.6)


def draw(shape):
    """TrueType glyph of a (multi)polygon: outer contours clockwise."""
    pen = TTGlyphPen(None)
    polygons = getattr(shape, "geoms", [shape])
    for polygon in polygons:
        polygon = orient(polygon, sign=-1.0)
        for ring in [polygon.exterior, *polygon.interiors]:
            points = [(round(x), round(y)) for x, y in ring.coords[:-1]]
            deduped = [p for i, p in enumerate(points)
                       if p != points[i - 1]]
            if len(deduped) < 3:
                continue
            pen.moveTo(deduped[0])
            for p in deduped[1:]:
                pen.lineTo(p)
            pen.closePath()
    return pen.glyph()


def build(path):
    order = [".notdef", "space"]
    cmap = {ord(" "): "space"}
    glyphs = {}
    metrics = {}

    box = TTGlyphPen(None)
    for ring in ([(60, 0), (60, 700), (535, 700), (535, 0)],
                 [(110, 50), (485, 50), (485, 650), (110, 650)]):
        box.moveTo(ring[0])
        for p in ring[1:]:
            box.lineTo(p)
        box.closePath()
    glyphs[".notdef"] = box.glyph()
    metrics[".notdef"] = (ADVANCE, 60)
    glyphs["space"] = TTGlyphPen(None).glyph()
    metrics["space"] = (ADVANCE, 0)

    for character, outline in GLYPHS.items():
        name = "uni%04X" % ord(character)
        shape = outline_shape(outline)
        order.append(name)
        cmap[ord(character)] = name
        glyphs[name] = draw(shape)
        metrics[name] = (ADVANCE, round(shape.bounds[0]))

    builder = FontBuilder(UNITS_PER_EM, isTTF=True)
    builder.setupGlyphOrder(order)
    builder.setupCharacterMap(cmap)
    builder.setupGlyf(glyphs)
    builder.setupHorizontalMetrics(metrics)
    builder.setupHorizontalHeader(ascent=ASCENT, descent=-DESCENT)
    builder.setupNameTable({
        "familyName": "FauxfoilPlot",
        "styleName": "Regular",
        "uniqueFontIdentifier": "FauxfoilPlot-Regular",
        "fullName": "FauxfoilPlot Regular",
        "psName": "FauxfoilPlot-Regular",
        "version": "Version 1.000",
        "copyright": "Copyright 2024-2026 Christopher A. Lupp",
        "licenseDescription": "Licensed under the Apache License, Version 2.0",
        "licenseInfoURL": "http://www.apache.org/licenses/LICENSE-2.0",
    })
    builder.setupOS2(sTypoAscender=ASCENT, sTypoDescender=-DESCENT,
                     sTypoLineGap=0, usWinAscent=ASCENT, usWinDescent=DESCENT,
                     sCapHeight=CAP_HEIGHT, sxHeight=round(4 * GRID),
                     fsType=0, achVendID="SCIR")
    builder.setupPost(isFixedPitch=1)
    builder.save(path)


if __name__ == "__main__":
    build(pathlib.Path(__file__).with_name("FauxfoilPlot.ttf"))
