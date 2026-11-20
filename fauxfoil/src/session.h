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

/// @file session.h
/// The fauXFOIL command interpreter: XFOIL's menus, prompts and output, with
/// Scirocco doing the analysis.

#include <iosfwd>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "airfoil.h"
#include "display.h"
#include "plots.h"
#include "polar.h"
#include "scirocco/config.h"
#include "scirocco/data/param.h"

namespace scirocco {
class Scirocco;
struct FoilState;
}  // namespace scirocco

namespace fauxfoil {

/// One interactive session reading commands from @p in and writing to
/// @p out, as XFOIL does on a terminal or from a piped script.
class Session {
   public:
    /// @param version Version string shown in the banner and polar files.
    Session(std::istream& in, std::ostream& out, std::string version);
    ~Session();
    Session(const Session&) = delete;
    Session& operator=(const Session&) = delete;

    /// Run the top-level menu until QUIT or end of input.
    /// @param args Command-line arguments; a first argument is loaded as an
    ///        airfoil file, as "xfoil file.dat" does.
    /// @return Process exit status.
    int Run(const std::vector<std::string>& args);

    /// Show plots in @p display (not owned); null leaves the session
    /// text-only. Plots are still built for hardcopies either way.
    void SetDisplay(Display* display);

    /// Text printed after the banner, such as why plots are unavailable.
    void SetNotice(std::string notice) { notice_ = std::move(notice); }

   private:
    enum class Target { Alpha, Cl, InviscidCl };

    /* Prompts. Each throws on end of input. */
    void Prompt(const std::string& label, const std::string& kind);
    std::string ReadLine();
    std::string AskCommand(const std::string& label);
    void AskReal(const std::string& prompt, double& value);
    void AskInt(const std::string& prompt, int& value);
    std::string AskString(const std::string& prompt);
    bool AskYesNo(const std::string& prompt);
    bool ConfirmOverwrite(const std::string& path, const char* refused);

    /* Menus */
    void TopLevel();
    void Oper();
    void Vpar();
    void Gdes();
    void Ppar();
    void Plop();
    [[nodiscard]] std::string OperLabel() const;

    /* Geometry */
    void Naca(const std::string& arg);
    void Load(const std::string& arg);
    void Save(const std::string& arg, bool labeled);
    void Pane();
    void SetCurrentFromBuffer();
    void PrintThicknessCamber(const MatrixXd& x);
    void PrintTrailingEdge();
    void PrintPanelingUsed();
    void Rename(const std::string& arg, bool increment);
    [[nodiscard]] bool BufferMatchesCurrent() const;

    /* Analysis */
    bool PolarLocked();
    void ApplyParameters(scirocco::FoilState& state) const;
    void SolvePoint(Target target, double value);
    void Sequence(Target target, const std::vector<double>& numbers);
    void OnIteration(const scirocco::FoilState& state, int iteration);
    void PrintMachRe() const;
    void PrintSonic(bool extra_blank) const;
    void ShowVpar() const;

    /* Polars and output files */
    void TogglePacc();
    void AddPolarPoint();
    void WritePolarFile(const Polar& polar);
    void ListPolars(int index) const;
    void PolarSummary() const;
    void WritePolar(const std::vector<double>& numbers);
    void DeletePolar(const std::vector<double>& numbers);
    void ReadPolar(const std::string& arg);
    void WriteCp(const std::string& arg);
    void WriteDump(const std::string& arg);
    void MinimumCp() const;

    /* Plots */
    enum class PlotKind {
        None,
        Cp,
        CpSequence,
        CpVectors,
        Bl,
        Polar,
        PolarAirfoils,
        Airfoil,
        Curvature
    };
    void Show(Scene scene, PlotKind kind);
    [[nodiscard]] bool HaveSolution() const;
    [[nodiscard]] Solution CurrentSolution() const;
    void PlotCp();
    void PlotPoint();
    void PlotCpVectors();
    void PlotBl();
    void PlotPolars();
    void PlotAirfoil();
    void Replot();
    void Hardcopy();
    void Vplo();
    void AskAxis(const char* name, Axis& axis,
                 const std::vector<double>& given);

    /* Cursor-driven and reference-data plot commands */
    bool CursorReady();
    std::optional<CursorHit> Cursor(bool done_button);
    bool MarkCorners(Point& a, Point& b, bool absolute);
    void Reshow();
    bool ZoomCommand(const std::string& name);
    void SizeCommand(double& size, const std::vector<double>& given);
    void Annotate();
    void Profiles(const std::vector<double>& given);
    void ProfilesAtCursor();
    void ProfileWeight(const std::vector<double>& given);
    void ReadCpReference();
    void ReadForceReference();
    void ReadReferencePolar(const std::string& arg);
    void DeleteReferencePolar(const std::vector<double>& given);
    void PolarAxesCommand();
    void BlBlowup();
    void DumpBlVariable(const std::string& arg);
    void OverlayBlVariable(const std::string& arg);
    void GdesBlowup();
    void GdesWindow();
    void GdesOverlay(const std::string& arg);
    void Curvatures();
    void Distance();
    void TickSize(const std::vector<double>& given);

    std::istream& in_;
    std::ostream& out_;
    std::string version_;
    std::string notice_;

    /* Buffer airfoil, in Scirocco's clockwise order */
    MatrixXd buffer_;
    std::string name_ = "";
    bool normalize_ = false;
    bool write_clockwise_ = false;
    scirocco::Vector2d xref_{0.25, 0.0};

    /* Current (paneled) airfoil */
    std::unique_ptr<scirocco::Scirocco> foil_;
    int nodes_ = 160;
    bool announce_vorticity_ = true;

    /* Operating point and boundary-layer settings */
    bool viscous_ = false;
    double re_ = 0.0;
    double mach_ = 0.0;
    bool bl_initialized_ = false;
    scirocco::Param bl_;
    double wake_length_ = 1.0;
    double vacc_ = 0.01;
    std::string last_solve_;
    int iteration_ = 0;

    /* Polars */
    std::vector<Polar> polars_;
    bool pacc_ = false;
    int active_ = -1;
    bool graphics_ = false;

    /* Plots */
    Display* display_ = nullptr;
    Scene plot_;
    PlotKind plot_kind_ = PlotKind::None;
    std::vector<Scene> hardcopies_;
    CpOptions cp_options_;
    bool sequence_polar_ = false;
    bool in_sequence_ = false;
    std::vector<Solution> sequence_;
    BlVariable bl_variable_ = BlVariable::Hk;
    BlAxes bl_axes_;
    AirfoilOptions airfoil_options_;

    /* Zoom, plot-object sizes and annotation settings */
    Rect zoom_;
    double oper_size_ = 10.0;
    double vplo_size_ = 10.0;
    double gdes_size_ = 10.0;
    double anno_width_ = 1.0;
    int anno_symbol_ = 0;
    double profile_weight_ = 0.02;

    /* Reference data */
    bool cp_reference_on_ = false;
    std::string cp_reference_file_;
    std::optional<CpReference> cp_reference_;
    bool force_reference_on_ = false;
    std::optional<ForceReference> force_reference_;
    std::vector<ReferencePolar> reference_polars_;
    PolarAxes polar_axes_;
};

}  // namespace fauxfoil
