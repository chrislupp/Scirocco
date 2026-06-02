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

    Philote-Cpp MDO Discipline Interface
*/
#pragma once

#include <explicit.h>

#include <memory>
#include <string>

#include "scirocco/scirocco.h"

#ifdef SCIROCCO_HDF5
#include "scirocco/io/hdf5_writer.h"
#endif

namespace scirocco {

/**
 * @brief Philote-Cpp explicit discipline wrapping the Scirocco airfoil solver.
 *
 * Exposes inviscid/viscous airfoil analysis as an MDO discipline with
 * analytical adjoint gradients for panel node coordinates and operating
 * conditions.
 *
 * The airfoil topology (NACA digits or initial coordinates) is defined via
 * options. The optimizer controls panel node positions and operating conditions
 * as continuous inputs.
 *
 * The discipline lifecycle follows Philote's convention:
 *   Initialize -> SetOptions -> Configure -> Setup -> SetupPartials ->
 *   Compute/ComputePartials (repeated)
 *
 * @par Options
 * - `naca` (string): NACA digit string (e.g., "2412")
 * - `coordinates` (dict): Airfoil coordinates {x: [...], z: [...]}
 * - `npanel` (int): Number of panels (default 199)
 *
 * One of `naca` or `coordinates` must be provided.
 *
 * @par Inputs
 * - `x` ({2*N}): Panel node coordinates, flattened row-major from [2 x N]
 * - `alpha` ({1}): Angle of attack [deg]
 * - `Re` ({1}): Reynolds number
 * - `Ma` ({1}): Mach number
 *
 * @par Outputs
 * - `cl` ({1}): Lift coefficient
 * - `cd` ({1}): Drag coefficient
 * - `cm` ({1}): Moment coefficient
 */
class SiroccoDiscipline : public philote::ExplicitDiscipline {
   public:
    SiroccoDiscipline() = default;
    ~SiroccoDiscipline() noexcept = default;

    void Initialize() override;
    void SetOptions(const google::protobuf::Struct& options) override;
    void Configure() override;
    void Setup() override;
    void SetupPartials() override;
    void Compute(const philote::Variables& inputs,
                 philote::Variables& outputs) override;
    void ComputePartials(const philote::Variables& inputs,
                         philote::Partials& partials) override;

   private:
    std::unique_ptr<Scirocco> solver_;
    int npanel_ = 199;
    int nnodes_ = 0;

    // ComputeSensitivities can throw on rare intermediate geometries even
    // with PTC stabilizing the primal (a singular adjoint Jacobian R_V is
    // independent of primal convergence). To keep the optimizer making
    // progress we cache the last successful Sensitivity here and reuse it
    // on failure. The first failure with no cache returns zeros and
    // prints a warning.
    Sensitivity cached_sens_;
    bool has_cached_sens_ = false;

#ifdef SCIROCCO_HDF5
    std::unique_ptr<Hdf5Writer> hdf5_writer_;
#endif
};

}  // namespace scirocco
