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

    Scirocco Philote Discipline Implementation
*/

#include "scirocco/philote_discipline.h"

#include <Eigen/Core>
#include <stdexcept>
#include <string>
#include <vector>

#include "option_validation.h"
#include "scirocco/core/geometry.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/core/paneling.h"

namespace scirocco {

void SiroccoDiscipline::Initialize() {
    AddOption("naca", "string");
    AddOption("coordinates", "dict");
    AddOption("npanel", "int");
    AddOption("hdf5_file", "string");
}

void SiroccoDiscipline::SetOptions(const google::protobuf::Struct& options) {
    solver_.reset();
    nnodes_ = 0;
    const auto& fields = options.fields();
    if ((fields.count("naca") && !fields.at("naca").has_string_value()) ||
        (fields.count("coordinates") &&
         !fields.at("coordinates").has_struct_value()) ||
        (fields.count("hdf5_file") &&
         !fields.at("hdf5_file").has_string_value()))
        throw std::invalid_argument(
            "Invalid type for naca, coordinates, or hdf5_file");

    if (fields.count("npanel")) {
        npanel_ = IntegerOption(fields.at("npanel"), "npanel", 4, kMaxPanels);
    }

    bool has_naca =
        fields.count("naca") && fields.at("naca").has_string_value();
    bool has_coords = fields.count("coordinates") &&
                      fields.at("coordinates").has_struct_value();

    if (has_naca && has_coords) {
        throw std::invalid_argument(
            "Only one of 'naca' or 'coordinates' may be specified");
    }

    if (!has_naca && !has_coords) {
        throw std::invalid_argument(
            "One of 'naca' or 'coordinates' must be specified");
    }

    if (has_naca) {
        std::string naca = fields.at("naca").string_value();
        solver_ = std::make_unique<Scirocco>(naca, npanel_);
    } else {
        const auto& coord_fields =
            fields.at("coordinates").struct_value().fields();
        if (!coord_fields.count("x") || !coord_fields.count("z") ||
            !coord_fields.at("x").has_list_value() ||
            !coord_fields.at("z").has_list_value())
            throw std::invalid_argument("coordinates requires x and z lists");
        const auto& x_list = coord_fields.at("x").list_value();
        const auto& z_list = coord_fields.at("z").list_value();

        int npts = x_list.values_size();
        if (npts < 4 || npts != z_list.values_size())
            throw std::invalid_argument("x and z must have equal lengths >= 4");
        MatrixXd coords(2, npts);
        for (int i = 0; i < npts; ++i) {
            coords(0, i) = NumberOption(x_list.values(i), "x coordinate");
            coords(1, i) = NumberOption(z_list.values(i), "z coordinate");
        }

        solver_ = std::make_unique<Scirocco>(coords, npanel_);
    }

    nnodes_ = static_cast<int>(solver_->GetFoilCoords().cols());

#ifdef SCIROCCO_HDF5
    if (fields.count("hdf5_file") &&
        fields.at("hdf5_file").has_string_value()) {
        hdf5_writer_ =
            std::make_unique<Hdf5Writer>(fields.at("hdf5_file").string_value());
    }
#endif
}

void SiroccoDiscipline::Configure() {}

void SiroccoDiscipline::Setup() {
    if (!solver_)
        throw std::runtime_error("Valid options are required before Setup");
    AddInput("x", {2 * nnodes_}, "m");
    AddInput("alpha", {1}, "deg");
    AddInput("Re", {1}, "");
    AddInput("Ma", {1}, "");

    AddOutput("cl", {1}, "");
    AddOutput("cd", {1}, "");
    AddOutput("cm", {1}, "");
}

void SiroccoDiscipline::SetupPartials() {
    std::vector<std::string> outputs = {"cl", "cd", "cm"};
    std::vector<std::string> inputs = {"x", "alpha", "Re", "Ma"};

    for (const auto& out : outputs) {
        for (const auto& in : inputs) {
            DeclarePartials(out, in);
        }
    }
}

void SiroccoDiscipline::Compute(const philote::Variables& inputs,
                                philote::Variables& outputs) {
    if (!solver_)
        throw std::runtime_error("Valid options are required before Compute");
    const auto& x_var = inputs.at("x");
    auto& state = solver_->GetState();
    // Invalidate the previous solution before any geometry update can fail.
    ClearSolution(state);
    if (x_var.Size() != 2 * static_cast<size_t>(nnodes_) ||
        inputs.at("alpha").Size() != 1 || inputs.at("Re").Size() != 1 ||
        inputs.at("Ma").Size() != 1)
        throw std::invalid_argument(
            "Input dimensions do not match the configured discipline");
    auto& foil_x = state.foil.x;
    for (int i = 0; i < nnodes_; ++i) {
        foil_x(0, i) = x_var(i);
        foil_x(1, i) = x_var(nnodes_ + i);
    }

    /* After overwriting panel node coordinates in place, refresh the
       derived geometric state (arc length, tangents) so it's consistent
       with the new foil.x. Without this the BL solver runs on the new
       nodes with stale tangents from the initial paneling, producing an
       inconsistent state and breaking the adjoint LU factorization on R_V. */
    RefreshPanelGeometry(state.foil);

    double alpha = inputs.at("alpha")(0);
    double Re = inputs.at("Re")(0);
    double Ma = inputs.at("Ma")(0);

    solver_->SetOper(alpha, -1.0, Re, Ma);
    if (!solver_->Solve()) {
        throw std::runtime_error("Scirocco analysis did not converge");
    }

    outputs.at("cl")(0) = solver_->GetCl();
    outputs.at("cd")(0) = solver_->GetCd();
    outputs.at("cm")(0) = solver_->GetCm();

#ifdef SCIROCCO_HDF5
    if (hdf5_writer_) {
        hdf5_writer_->WriteIteration(*solver_);
    }
#endif
}

void SiroccoDiscipline::ComputePartials(const philote::Variables& inputs,
                                        philote::Partials& partials) {
    if (!solver_)
        throw std::runtime_error(
            "Valid options are required before ComputePartials");
    const auto& state = solver_->GetState();
    bool same_inputs = inputs.at("alpha")(0) == state.oper.alpha &&
                       (state.oper.viscous ? inputs.at("Re")(0) == state.oper.Re
                                           : inputs.at("Re")(0) <= 0.0) &&
                       inputs.at("Ma")(0) == state.oper.Ma;
    for (int i = 0; i < nnodes_; ++i) {
        same_inputs = same_inputs && inputs.at("x")(i) == state.foil.x(0, i) &&
                      inputs.at("x")(nnodes_ + i) == state.foil.x(1, i);
    }
    if (!same_inputs) {
        throw std::runtime_error(
            "ComputePartials requires Compute at the same inputs");
    }
    const Sensitivity sens = solver_->ComputeSensitivities();

    /* d(cl,cd,cm)/dx — flatten [2 x N] row-major into [2N] */
    auto& dcl_dx = partials[{"cl", "x"}];
    auto& dcd_dx = partials[{"cd", "x"}];
    auto& dcm_dx = partials[{"cm", "x"}];

    for (int i = 0; i < nnodes_; ++i) {
        dcl_dx(i) = sens.dcl_dx(0, i);
        dcl_dx(nnodes_ + i) = sens.dcl_dx(1, i);

        dcd_dx(i) = sens.dcd_dx(0, i);
        dcd_dx(nnodes_ + i) = sens.dcd_dx(1, i);

        dcm_dx(i) = sens.dcm_dx(0, i);
        dcm_dx(nnodes_ + i) = sens.dcm_dx(1, i);
    }

    /* d(cl,cd,cm)/d(alpha) */
    partials[{"cl", "alpha"}](0) = sens.dcl_dalpha;
    partials[{"cd", "alpha"}](0) = sens.dcd_dalpha;
    partials[{"cm", "alpha"}](0) = sens.dcm_dalpha;

    /* d(cl,cd,cm)/d(Re) */
    partials[{"cl", "Re"}](0) = sens.dcl_dRe;
    partials[{"cd", "Re"}](0) = sens.dcd_dRe;
    partials[{"cm", "Re"}](0) = sens.dcm_dRe;

    /* d(cl,cd,cm)/d(Ma) */
    partials[{"cl", "Ma"}](0) = sens.dcl_dMa;
    partials[{"cd", "Ma"}](0) = sens.dcd_dMa;
    partials[{"cm", "Ma"}](0) = sens.dcm_dMa;
}

}  // namespace scirocco
