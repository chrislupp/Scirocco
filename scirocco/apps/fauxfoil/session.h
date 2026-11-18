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
/// The fauXfoil command interpreter: XFOIL's menus, prompts and output, with
/// Scirocco doing the analysis.

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

#include "airfoil.h"
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

    std::istream& in_;
    std::ostream& out_;
    std::string version_;

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
};

}  // namespace fauxfoil
