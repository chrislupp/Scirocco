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
    int node = -1;  ///< Airfoil node, or -1 in the wake
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

/// Measured Cp(x) overlaid on the Cp plot (CREF).
struct CpReference {
    std::string label;
    std::vector<Point> points;  ///< (x, Cp)
};

/// Measured forces listed beside the computed ones (FREF).
struct ForceReference {
    std::string label;
    double mach = 0.0;
    double re = 0.0;
    double alpha = 0.0;
    double cl = 0.0;
    double cd = 0.0;
    double cm = 0.0;
};

/// Settings shared by the Cp(x) plots.
struct CpOptions {
    double cp_min = -2.0;  ///< Top of the Cp axis (CPMI)
    bool grid = false;     ///< Grid overlay (GRID)
    bool inviscid = true;  ///< Inviscid Cp overlay on viscous Cp (CPIO)
    std::string version;   ///< Program version for the corner tag
    const CpReference* cp_reference = nullptr;        ///< CREF data
    const ForceReference* force_reference = nullptr;  ///< FREF data
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
    bool vs_s = false;     ///< Plot against arc length instead of x (XPLO)
    bool grid = true;      ///< Grid (GRID)
    bool symbols = false;  ///< Node symbols on the curves (SYMB)
    bool labels = true;    ///< Header and axis titles (LABE)
    bool clip = false;     ///< Clip curves to the axes (CLIP)
    std::vector<Point> overlay;  ///< Data from a file (OVER), as (x, value)
};

/// Default axes a boundary-layer plot would use.
BlAxes DefaultBlAxes(const Solution& solution, BlVariable variable,
                     const BlAxes& current);

/// The curves a boundary-layer plot draws, as (x or s, value) points.
std::vector<std::vector<Point>> BlData(const Solution& solution,
                                       BlVariable variable, bool vs_s);

/// Axis over [lo, hi] with round limits and at most @p intervals labels.
Axis NiceAxisFor(double lo, double hi, int intervals);

/// Plot one boundary-layer variable on both surfaces and the wake.
Scene BlPlot(const Solution& solution, BlVariable variable, const BlAxes& axes);

/// Measured polar data overlaid on the polar plot (RGET).
struct ReferencePolar {
    std::string label;
    std::vector<Point> cd_cl;     ///< (CD, CL)
    std::vector<Point> alpha_cl;  ///< (alpha, CL)
    std::vector<Point> alpha_cm;  ///< (alpha, CM)
    std::vector<Point> xtr_cl;    ///< (xtr/c, CL)
};

/// Polar plot axis limits set with PPAX; unset axes cover the data.
struct PolarAxes {
    bool set = false;
    Axis alpha{-4.0, 10.0, 2.0};
    Axis cl{0.0, 1.5, 0.5};
    Axis cd{0.0, 0.04, 0.01};
    Axis cm{-0.25, 0.0, 0.05};  ///< As XFOIL lists it for "-CM"
};

/// The polar plot: CL-CD, CL-alpha with CM-alpha, and CL-xtr (PPLO).
Scene PolarPlot(const std::vector<Polar>& polars, const std::string& version,
                const PolarAxes& axes = {},
                const std::vector<ReferencePolar>& references = {});

/// Axes the polar plot uses for these polars.
PolarAxes DefaultPolarAxes(const std::vector<Polar>& polars,
                           const PolarAxes& axes);

/// The airfoils stored with polars, overlaid in the polars' colours
/// (APLO). @p which lists 1-based polar numbers.
Scene PolarAirfoilsPlot(const std::vector<Polar>& polars,
                        const std::vector<int>& which);

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
    bool grid = true;           ///< GRID
    bool ticks = false;         ///< Node tick marks (TICK)
    double tick_size = 0.0005;  ///< Tick length / perimeter (TSIZ)
    bool parameters = true;     ///< Geometric parameters (GPAR)
    bool inches = false;        ///< Unscaled, one inch per unit (INPL)
    bool window_set = false;    ///< Plot only this region (BLOW, WIND)
    Rect window;                ///< In airfoil coordinates
    /// Another airfoil drawn over the buffer airfoil (OVER).
    std::vector<Point> overlay;
    GeometryInfo overlay_info;
};

/// The buffer airfoil on a grid, as GDES shows it.
Scene AirfoilPlot(const std::vector<Point>& contour, const GeometryInfo& info,
                  const AirfoilOptions& options);

/// Data region the GDES plot shows by default.
Rect AirfoilWindow(const std::vector<Point>& contour);

/// Curvature of a contour at each point, from the circle through each
/// point and its neighbours.
std::vector<double> Curvature(const std::vector<Point>& contour);

/// Curvature (raised to @p exponent) against x on both surfaces (CPLO).
Scene CurvaturePlot(const std::vector<Point>& contour, double exponent);

/// Format a number with a fixed count of decimals, never as "-0.0".
std::string Fixed(double value, int decimals);

}  // namespace fauxfoil
