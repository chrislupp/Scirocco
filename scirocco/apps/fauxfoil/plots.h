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

/// @file plots.h
/// XFOIL's plots, built as scenes: Cp(x), Cp vectors, boundary-layer
/// variables, polars and the airfoil geometry.

#include <string>
#include <vector>

#include "plot_scene.h"
#include "polar.h"

namespace fauxfoil {

/// Boundary-layer values at one station.
struct BlStation {
    double x = 0.0;     ///< x/c
    double s = 0.0;     ///< Arc length from the stagnation point
    double hk = 0.0;    ///< Kinematic shape parameter
    double ds = 0.0;    ///< Displacement thickness
    double th = 0.0;    ///< Momentum thickness
    double ue = 0.0;    ///< Edge velocity / Vinf
    double cf = 0.0;    ///< Skin friction coefficient
    double cd = 0.0;    ///< Dissipation coefficient
    double sa = 0.0;    ///< Amplification (laminar) or sqrt(ctau) (turbulent)
    double cteq = 0.0;  ///< Equilibrium sqrt(ctau)
    double ret = 0.0;   ///< Re_theta
    bool turbulent = false;
};

/// What the plots need from one solved operating point.
struct Solution {
    std::string name;
    bool viscous = false;
    double alpha = 0.0;  ///< Degrees
    double cl = 0.0;
    double cm = 0.0;
    double cd = 0.0;
    double cdp = 0.0;  ///< Inviscid pressure drag (near field)
    double re = 0.0;
    double mach = 0.0;
    double ncrit = 9.0;
    double xtr_top = 1.0;
    double xtr_bot = 1.0;

    std::vector<Point> airfoil;  ///< Panel nodes, Scirocco (clockwise) order
    std::vector<double> cp;      ///< Cp at the nodes
    std::vector<double> cpi;     ///< Inviscid Cp at the nodes (viscous only)
    std::vector<int> upper;      ///< Node indices, stagnation to upper TE
    std::vector<int> lower;      ///< Node indices, stagnation to lower TE
    std::vector<Point> wake;     ///< Wake points (viscous only)
    std::vector<double> cp_wake;
    std::vector<double> cpi_wake;

    std::vector<BlStation> bl_top;     ///< Stagnation to upper TE
    std::vector<BlStation> bl_bottom;  ///< Stagnation to lower TE
    std::vector<BlStation> bl_wake;
    std::vector<Point> dstar_top;     ///< Displacement surface, upper side
    std::vector<Point> dstar_bottom;  ///< Displacement surface, lower side
    std::vector<Point> wake_upper;    ///< Wake displacement edges
    std::vector<Point> wake_lower;
};

/// Settings shared by the Cp(x) plots.
struct CpOptions {
    double cp_min = -2.0;  ///< Top of the Cp axis (CPMI)
    bool grid = false;     ///< Grid overlay (GRID)
    bool inviscid = true;  ///< Inviscid Cp overlay on viscous Cp (CPIO)
    std::string version;   ///< Program version for the corner tag
};

/// Cp(x) with the airfoil below and the operating point listed (CPX).
Scene CpPlot(const Solution& solution, const CpOptions& options);

/// Cp(x) of every point of a sequence, with a table of results.
Scene CpSequencePlot(const std::vector<Solution>& solutions,
                     const CpOptions& options);

/// Airfoil at its angle of attack with surface pressure vectors (CPV).
Scene CpVectorPlot(const Solution& solution);

/// Boundary-layer plot variables (VPLO).
enum class BlVariable {
    Hk,
    DstarThetaTop,
    DstarThetaBottom,
    Ue,
    Cf,
    Dissipation,
    Amplification,
    Ctau,
    RetTheta,
    LogRetTheta
};

/// An axis: limits and label spacing.
struct Axis {
    double min = 0.0;
    double max = 1.0;
    double del = 0.2;
};

/// Axis limits of a boundary-layer plot; unset limits are chosen from the
/// data.
struct BlAxes {
    bool x_set = false;
    Axis x{0.0, 1.4, 0.2};
    bool y_set = false;
    Axis y;
    bool vs_s = false;  ///< Plot against arc length instead of x (XPLO)
    bool grid = true;   ///< Grid (GRID)
};

/// Default axes a boundary-layer plot would use.
BlAxes DefaultBlAxes(const Solution& solution, BlVariable variable,
                     const BlAxes& current);

/// Plot one boundary-layer variable on both surfaces and the wake.
Scene BlPlot(const Solution& solution, BlVariable variable, const BlAxes& axes);

/// The polar plot: CL-CD, CL-alpha with CM-alpha, and CL-xtr (PPLO).
Scene PolarPlot(const std::vector<Polar>& polars, const std::string& version);

/// Geometric parameters listed on the GDES plot.
struct GeometryInfo {
    std::string name;
    double area = 0.0;
    double thickness = 0.0;
    double camber = 0.0;
    double le_radius = 0.0;
    double te_angle = 0.0;  ///< Degrees
};

/// Options of the GDES airfoil plot.
struct AirfoilOptions {
    bool grid = true;        ///< GRID
    bool ticks = false;      ///< Node tick marks (TICK)
    bool parameters = true;  ///< Geometric parameters (GPAR)
};

/// The buffer airfoil on a grid, as GDES shows it.
Scene AirfoilPlot(const std::vector<Point>& contour, const GeometryInfo& info,
                  const AirfoilOptions& options);

/// Format a number with a fixed count of decimals, never as "-0.0".
std::string Fixed(double value, int decimals);

}  // namespace fauxfoil
