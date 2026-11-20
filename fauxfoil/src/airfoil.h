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

/// @file airfoil.h
/// Airfoil coordinate files and the geometry summaries XFOIL prints.
///
/// Coordinates are [2 x N] matrices. Scirocco orders points clockwise (lower
/// trailing edge, leading edge, upper trailing edge); XFOIL files are
/// counterclockwise (upper trailing edge first).

#include <string>

#include "scirocco/config.h"

namespace fauxfoil {

using scirocco::MatrixXd;

/// Maximum thickness and camber and their chordwise locations.
struct ThicknessCamber {
    double thickness = 0.0;
    double x_thickness = 0.0;
    double camber = 0.0;
    double x_camber = 0.0;
};

/// Thickness and camber from the upper and lower surfaces, split at the
/// leading edge and measured in the chord-line frame.
ThicknessCamber GetThicknessCamber(const MatrixXd& x);

/// Index of the point of minimum x.
int LeadingEdgeIndex(const MatrixXd& x);

/// Leading edge as XFOIL's LEFIND defines it: the point on a cubic spline
/// through the points where the surface is normal to the line from the
/// trailing-edge midpoint.
scirocco::Vector2d LeadingEdge(const MatrixXd& x);

/// Distance between the first and last points.
double TrailingEdgeGap(const MatrixXd& x);

/// Chord: distance from the leading edge to the trailing-edge midpoint.
double Chord(const MatrixXd& x);

/// True if the points run counterclockwise (positive enclosed area).
bool IsCounterclockwise(const MatrixXd& x);

/// The points in reverse order.
MatrixXd Reversed(const MatrixXd& x);

/// Translate and scale to unit chord with the leading edge at the origin.
MatrixXd Normalized(const MatrixXd& x);

/// Largest angle [deg] between adjacent panels, and the node where it occurs.
double MaxCornerAngle(const MatrixXd& x, int* node);

/// Contents of an airfoil coordinate file.
struct AirfoilFile {
    bool labeled = false; /**< First line holds the name */
    std::string name;
    MatrixXd x; /**< Points in file order */
};

/// Read a labeled, plain, ISES or MSES single-element coordinate file.
/// @return false with @p error set if the file cannot be read.
bool ReadAirfoilFile(const std::string& path, AirfoilFile& file,
                     std::string& error);

/// Write points in the given order, one "x y" pair per line in XFOIL's
/// (1X,2G15.7) layout, after a name line if @p name is non-null.
bool WriteAirfoilFile(const std::string& path, const std::string* name,
                      const MatrixXd& x);

}  // namespace fauxfoil
