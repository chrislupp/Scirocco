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
#include "scirocco/io/hdf5_writer.h"

#include <H5Cpp.h>

#include <chrono>
#include <cstdio>

#include "scirocco/data/sensitivity.h"
#include "scirocco/data/state.h"
#include "scirocco/scirocco.h"

namespace scirocco {

namespace {

void WriteScalar(H5::Group& grp, const std::string& name, double value) {
    H5::DataSpace space(H5S_SCALAR);
    auto ds = grp.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space);
    ds.write(&value, H5::PredType::NATIVE_DOUBLE);
}

void WriteInt(H5::Group& grp, const std::string& name, int value) {
    H5::DataSpace space(H5S_SCALAR);
    auto ds = grp.createDataSet(name, H5::PredType::NATIVE_INT, space);
    ds.write(&value, H5::PredType::NATIVE_INT);
}

void WriteVector(H5::Group& grp, const std::string& name, const VectorXd& vec) {
    if (vec.size() == 0) return;
    hsize_t dims[1] = {static_cast<hsize_t>(vec.size())};
    H5::DataSpace space(1, dims);
    auto ds = grp.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space);
    ds.write(vec.data(), H5::PredType::NATIVE_DOUBLE);
}

void WriteMatrix(H5::Group& grp, const std::string& name, const MatrixXd& mat) {
    if (mat.size() == 0) return;
    /* Eigen is column-major; write as (rows x cols) and transpose to row-major
       for HDF5 compatibility with h5py/numpy (which expect row-major C order).
     */
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        row_major = mat;
    hsize_t dims[2] = {static_cast<hsize_t>(mat.rows()),
                       static_cast<hsize_t>(mat.cols())};
    H5::DataSpace space(2, dims);
    auto ds = grp.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space);
    ds.write(row_major.data(), H5::PredType::NATIVE_DOUBLE);
}

void WriteStringAttr(H5::Group& grp, const std::string& name,
                     const std::string& value) {
    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
    H5::DataSpace space(H5S_SCALAR);
    auto attr = grp.createAttribute(name, str_type, space);
    attr.write(str_type, value);
}

std::string CurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S",
                  std::gmtime(&time_t_now));
    return std::string(buf);
}

}  // namespace

struct Hdf5Writer::Impl {
    H5::H5File file;
    int iteration_count = 0;
    bool is_open = false;

    void WriteOperatingConditions(H5::Group& iter_grp, const FoilState& state) {
        auto grp = iter_grp.createGroup("operating_conditions");
        WriteScalar(grp, "alpha", state.oper.alpha);
        WriteScalar(grp, "Re", state.oper.Re);
        WriteScalar(grp, "Ma", state.oper.Ma);
        WriteScalar(grp, "Vinf", state.oper.Vinf);
        WriteScalar(grp, "rho", state.oper.rho);
        WriteInt(grp, "viscous", state.oper.viscous ? 1 : 0);
    }

    void WriteGeometry(H5::Group& iter_grp, const FoilState& state) {
        auto grp = iter_grp.createGroup("geometry");
        WriteScalar(grp, "chord", state.geom.chord);
        WriteScalar(grp, "hTE", state.geom.hTE);
        WriteMatrix(grp, "xpoint", state.geom.xpoint);

        /* xref as a small vector */
        VectorXd xref(2);
        xref(0) = state.geom.xref(0);
        xref(1) = state.geom.xref(1);
        WriteVector(grp, "xref", xref);

        WriteStringAttr(grp, "name", state.geom.name);
    }

    void WritePanels(H5::Group& iter_grp, const FoilState& state) {
        auto grp = iter_grp.createGroup("panels");

        /* Foil panels */
        auto foil_grp = grp.createGroup("foil");
        WriteInt(foil_grp, "N", state.foil.N);
        WriteMatrix(foil_grp, "x", state.foil.x);
        WriteVector(foil_grp, "s", state.foil.s);
        WriteMatrix(foil_grp, "t", state.foil.t);

        /* Wake panels */
        auto wake_grp = grp.createGroup("wake");
        WriteInt(wake_grp, "N", state.wake.N);
        WriteMatrix(wake_grp, "x", state.wake.x);
        WriteVector(wake_grp, "s", state.wake.s);
        WriteMatrix(wake_grp, "t", state.wake.t);
    }

    void WriteCoefficients(H5::Group& iter_grp, const FoilState& state) {
        auto grp = iter_grp.createGroup("coefficients");
        WriteScalar(grp, "cl", state.post.cl);
        WriteScalar(grp, "cd", state.post.cd);
        WriteScalar(grp, "cdf", state.post.cdf);
        WriteScalar(grp, "cdp", state.post.cdp);
        WriteScalar(grp, "cdpi", state.post.cdpi);
        WriteScalar(grp, "cm", state.post.cm);
        WriteInt(grp, "converged", state.glob.conv ? 1 : 0);
    }

    void WriteDistributions(H5::Group& iter_grp, const FoilState& state) {
        auto grp = iter_grp.createGroup("distributions");
        WriteVector(grp, "cp", state.post.cp);
        WriteVector(grp, "cpi", state.post.cpi);
        WriteVector(grp, "cf", state.post.cf);
        WriteVector(grp, "ue", state.post.ue);
        WriteVector(grp, "uei", state.post.uei);
        WriteVector(grp, "th", state.post.th);
        WriteVector(grp, "ds", state.post.ds);
        WriteVector(grp, "sa", state.post.sa);
        WriteVector(grp, "Ret", state.post.Ret);
        WriteVector(grp, "Hk", state.post.Hk);
    }

    void WriteSensitivities(H5::Group& iter_grp, const Sensitivity& sens) {
        auto grp = iter_grp.createGroup("sensitivities");

        /* Shape sensitivities */
        auto shape_grp = grp.createGroup("shape");
        WriteMatrix(shape_grp, "dcl_dx", sens.dcl_dx);
        WriteMatrix(shape_grp, "dcd_dx", sens.dcd_dx);
        WriteMatrix(shape_grp, "dcm_dx", sens.dcm_dx);
        WriteMatrix(shape_grp, "dcdpi_dx", sens.dcdpi_dx);

        /* Scalar sensitivities */
        auto scalar_grp = grp.createGroup("scalar");
        WriteScalar(scalar_grp, "dcl_dalpha", sens.dcl_dalpha);
        WriteScalar(scalar_grp, "dcd_dalpha", sens.dcd_dalpha);
        WriteScalar(scalar_grp, "dcm_dalpha", sens.dcm_dalpha);
        WriteScalar(scalar_grp, "dcl_dMa", sens.dcl_dMa);
        WriteScalar(scalar_grp, "dcd_dMa", sens.dcd_dMa);
        WriteScalar(scalar_grp, "dcm_dMa", sens.dcm_dMa);
        WriteScalar(scalar_grp, "dcl_dRe", sens.dcl_dRe);
        WriteScalar(scalar_grp, "dcd_dRe", sens.dcd_dRe);
        WriteScalar(scalar_grp, "dcm_dRe", sens.dcm_dRe);
    }

    void UpdateNumIterationsAttr() {
        if (!is_open) return;
        H5::DataSpace space(H5S_SCALAR);
        if (file.attrExists("num_iterations")) {
            auto attr = file.openAttribute("num_iterations");
            attr.write(H5::PredType::NATIVE_INT, &iteration_count);
        } else {
            auto attr = file.createAttribute("num_iterations",
                                             H5::PredType::NATIVE_INT, space);
            attr.write(H5::PredType::NATIVE_INT, &iteration_count);
        }
    }
};

Hdf5Writer::Hdf5Writer(const std::string& filepath) : pimpl_(new Impl) {
    pimpl_->file = H5::H5File(filepath, H5F_ACC_TRUNC);
    pimpl_->is_open = true;
    pimpl_->iteration_count = 0;

    /* Write root metadata */
    H5::DataSpace space(H5S_SCALAR);
    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);

    std::string version = "0.1.0";
    auto ver_attr =
        pimpl_->file.createAttribute("scirocco_version", str_type, space);
    ver_attr.write(str_type, version);

    std::string ts = CurrentTimestamp();
    auto ts_attr = pimpl_->file.createAttribute("created", str_type, space);
    ts_attr.write(str_type, ts);

    pimpl_->UpdateNumIterationsAttr();
}

Hdf5Writer::~Hdf5Writer() { Close(); }

void Hdf5Writer::WriteIteration(const Scirocco& solver,
                                const Sensitivity* sensitivity) {
    WriteIteration(solver.GetState(), sensitivity);
}

void Hdf5Writer::WriteIteration(const FoilState& state,
                                const Sensitivity* sensitivity) {
    char name[32];
    std::snprintf(name, sizeof(name), "iteration_%03d",
                  pimpl_->iteration_count);

    auto iter_grp = pimpl_->file.createGroup(name);

    pimpl_->WriteOperatingConditions(iter_grp, state);
    pimpl_->WriteGeometry(iter_grp, state);
    pimpl_->WritePanels(iter_grp, state);
    pimpl_->WriteCoefficients(iter_grp, state);
    pimpl_->WriteDistributions(iter_grp, state);

    if (sensitivity) {
        pimpl_->WriteSensitivities(iter_grp, *sensitivity);
    }

    pimpl_->iteration_count++;
    pimpl_->UpdateNumIterationsAttr();
    pimpl_->file.flush(H5F_SCOPE_LOCAL);
}

void Hdf5Writer::Close() {
    if (pimpl_ && pimpl_->is_open) {
        pimpl_->UpdateNumIterationsAttr();
        pimpl_->file.close();
        pimpl_->is_open = false;
    }
}

int Hdf5Writer::NumIterations() const { return pimpl_->iteration_count; }

}  // namespace scirocco
