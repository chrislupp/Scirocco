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

    Philote Discipline Tests

    Validates that the Philote discipline interface produces the same primal
    outputs and analytical partials as standalone Scirocco, and that the
    partials match the mfoil complex-step reference data.
*/

#include "scirocco/philote_discipline.h"

#include <H5Cpp.h>
#include <google/protobuf/struct.pb.h>
#include <gtest/gtest.h>
#include <variable.h>

#include <Eigen/Core>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "scirocco/kulfan_discipline.h"
#include "scirocco/scirocco.h"

using namespace Eigen;

/* Tolerances matching sensitivity_reference_test.cpp */
constexpr double TOL_INVISCID = 1e-10;
/* ---------------------------------------------------------------------------
   HDF5 helpers (same as sensitivity_reference_test.cpp)
   ---------------------------------------------------------------------------
 */
template <typename T>
T readScalar(H5::H5File& file, const std::string& path) {
    H5::DataSet dataset = file.openDataSet(path);
    T value;
    // Use the destination memory type, not the on-disk storage width.
    static_assert(std::is_same_v<T, int> || std::is_same_v<T, double>);
    if constexpr (std::is_same_v<T, int>) {
        dataset.read(&value, H5::PredType::NATIVE_INT);
    } else {
        dataset.read(&value, H5::PredType::NATIVE_DOUBLE);
    }
    return value;
}

MatrixXd readMatrix(H5::H5File& file, const std::string& path) {
    H5::DataSet dataset = file.openDataSet(path);
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims);
    std::vector<double> buffer(dims[0] * dims[1]);
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
    MatrixXd result(dims[0], dims[1]);
    for (hsize_t i = 0; i < dims[0]; ++i)
        for (hsize_t j = 0; j < dims[1]; ++j)
            result(i, j) = buffer[i * dims[1] + j];
    return result;
}

std::vector<int> readIntVector(H5::H5File& file, const std::string& path) {
    H5::DataSet dataset = file.openDataSet(path);
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims);
    std::vector<int64_t> buffer(dims[0]);
    dataset.read(buffer.data(), H5::PredType::NATIVE_INT64);
    std::vector<int> result(dims[0]);
    for (size_t i = 0; i < dims[0]; ++i)
        result[i] = static_cast<int>(buffer[i]);
    return result;
}

/* ---------------------------------------------------------------------------
   Protobuf helpers
   ---------------------------------------------------------------------------
 */
namespace {

google::protobuf::Struct BuildNacaOptions(const std::string& naca, int npanel) {
    google::protobuf::Struct opts;
    (*opts.mutable_fields())["naca"].set_string_value(naca);
    (*opts.mutable_fields())["npanel"].set_number_value(npanel);
    return opts;
}

google::protobuf::Struct BuildCoordOptions(const scirocco::MatrixXd& coords,
                                           int npanel) {
    google::protobuf::Struct opts;
    (*opts.mutable_fields())["npanel"].set_number_value(npanel);

    auto& coord_struct =
        *(*opts.mutable_fields())["coordinates"].mutable_struct_value();
    auto& x_list = *(*coord_struct.mutable_fields())["x"].mutable_list_value();
    auto& z_list = *(*coord_struct.mutable_fields())["z"].mutable_list_value();

    for (int i = 0; i < coords.cols(); ++i) {
        x_list.add_values()->set_number_value(coords(0, i));
        z_list.add_values()->set_number_value(coords(1, i));
    }

    return opts;
}

template <typename Discipline>
void RunLifecycle(Discipline& disc, const google::protobuf::Struct& opts,
                  philote::Variables& inputs, philote::Variables& outputs,
                  philote::Partials& partials) {
    disc.Initialize();
    disc.SetOptions(opts);
    disc.Setup();
    disc.SetupPartials();

    for (const auto& meta : disc.var_meta()) {
        std::vector<size_t> shape;
        for (auto dim : meta.shape()) shape.push_back(static_cast<size_t>(dim));

        if (meta.type() == philote::kInput)
            inputs[meta.name()] = philote::Variable(philote::kInput, shape);
        else if (meta.type() == philote::kOutput)
            outputs[meta.name()] = philote::Variable(philote::kOutput, shape);
    }

    for (const auto& meta : disc.partials_meta()) {
        std::vector<size_t> shape;
        for (auto dim : meta.shape()) shape.push_back(static_cast<size_t>(dim));

        partials[std::make_pair(meta.name(), meta.subname())] =
            philote::Variable(philote::kPartial, shape);
    }
}

void SetCoordinateInputs(philote::Variables& inputs,
                         const scirocco::MatrixXd& foil_x, double alpha,
                         double Re, double Ma) {
    int nnodes = static_cast<int>(foil_x.cols());
    auto& x_var = inputs.at("x");
    for (int i = 0; i < nnodes; ++i) {
        x_var(i) = foil_x(0, i);
        x_var(nnodes + i) = foil_x(1, i);
    }
    inputs.at("alpha")(0) = alpha;
    inputs.at("Re")(0) = Re;
    inputs.at("Ma")(0) = Ma;
}

} /* namespace */

/* ===========================================================================
   Primal output tests: discipline vs. standalone Scirocco
   ===========================================================================
 */

TEST(PhilotePrimalTest, InviscidMatchesStandalone) {
    scirocco::SiroccoDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildNacaOptions("0012", 199), inputs, outputs,
                 partials);

    scirocco::Scirocco ref("0012", 199);
    ref.SetOper(5.0);
    ref.Solve();

    SetCoordinateInputs(inputs, ref.GetFoilCoords(), 5.0, 0.0, 0.0);
    disc.Compute(inputs, outputs);

    EXPECT_NEAR(outputs.at("cl")(0), ref.GetCl(), 1e-12);
    EXPECT_NEAR(outputs.at("cd")(0), ref.GetCd(), 1e-12);
    EXPECT_NEAR(outputs.at("cm")(0), ref.GetCm(), 1e-12);
}

TEST(PhilotePrimalTest, ViscousMatchesStandalone) {
    scirocco::SiroccoDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildNacaOptions("0012", 199), inputs, outputs,
                 partials);

    scirocco::Scirocco ref("0012", 199);
    ref.SetOper(3.0, -1.0, 1e6, 0.0);
    ref.Solve();

    SetCoordinateInputs(inputs, ref.GetFoilCoords(), 3.0, 1e6, 0.0);
    disc.Compute(inputs, outputs);

    /* The discipline now refreshes panel.s/panel.t from the input coords
       (via a single Spline2DCreate pass) while the reference Scirocco's
       panel derived state came from the constructor's two-pass
       SplineCurvature pipeline. Both are internally consistent but differ
       at the spline-numerical level: cl/cd agree to ~1e-9, cm to ~1e-8. */
    EXPECT_NEAR(outputs.at("cl")(0), ref.GetCl(), 1e-8);
    EXPECT_NEAR(outputs.at("cd")(0), ref.GetCd(), 1e-8);
    EXPECT_NEAR(outputs.at("cm")(0), ref.GetCm(), 1e-8);
}

TEST(PhilotePrimalTest, TutorialReynoldsSweep) {
    // Match examples/tutorials/tutorial_airfoil.py, including its finite TE.
    scirocco::KulfanDiscipline geometry;
    google::protobuf::Struct options;
    for (const auto* name : {"n_cst_upper", "n_cst_lower"})
        (*options.mutable_fields())[name].set_number_value(6);
    (*options.mutable_fields())["npanel"].set_number_value(100);
    (*options.mutable_fields())["te_upper"].set_number_value(0.00063);
    (*options.mutable_fields())["te_lower"].set_number_value(-0.00063);
    philote::Variables geometry_inputs, geometry_outputs;
    philote::Partials geometry_partials;
    RunLifecycle(geometry, options, geometry_inputs, geometry_outputs,
                 geometry_partials);
    const std::array<double, 6> coefficients{
        0.1706635358053721,  0.15626369325824663, 0.1546327069833401,
        0.14462202838023833, 0.13318428664336585, 0.15537328876198608};
    for (int i = 0; i < 6; ++i) {
        geometry_inputs.at("cst_upper")(i) = coefficients[i];
        geometry_inputs.at("cst_lower")(i) = coefficients[i];
    }
    geometry.Compute(geometry_inputs, geometry_outputs);
    scirocco::SiroccoDiscipline analysis;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(analysis, BuildNacaOptions("0012", 100), inputs, outputs,
                 partials);
    for (int i = 0; i < 202; ++i)
        inputs.at("x")(i) = geometry_outputs.at("x")(i);
    inputs.at("alpha")(0) = 2.0;
    inputs.at("Ma")(0) = 0.0;
    // Rounded tutorial output: catch nonconvergence and stale cached forces.
    const std::array<std::array<double, 3>, 4> cases{{{2.5e5, 0.27840, 0.00967},
                                                      {5e5, 0.21612, 0.00698},
                                                      {1e6, 0.20489, 0.00565},
                                                      {2e6, 0.21067, 0.00531}}};
    for (const auto& row : cases) {
        SCOPED_TRACE(row[0]);
        inputs.at("Re")(0) = row[0];
        ASSERT_NO_THROW(analysis.Compute(inputs, outputs));
        EXPECT_NEAR(outputs.at("cl")(0), row[1], 1e-5);
        EXPECT_NEAR(outputs.at("cd")(0), row[2], 1e-5);
        EXPECT_TRUE(std::isfinite(outputs.at("cm")(0)));
    }
}

/* ===========================================================================
   Coordinate-input options path
   ===========================================================================
 */

TEST(PhilotePrimalTest, CoordinateOptionsMatchNaca) {
    scirocco::Scirocco naca_src("0012", 99);
    auto coords = naca_src.GetFoilCoords();

    scirocco::Scirocco ref(coords, 99);
    ref.SetOper(5.0);
    ref.Solve();

    scirocco::SiroccoDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildCoordOptions(coords, 99), inputs, outputs,
                 partials);

    SetCoordinateInputs(inputs, ref.GetFoilCoords(), 5.0, 0.0, 0.0);
    disc.Compute(inputs, outputs);

    EXPECT_NEAR(outputs.at("cl")(0), ref.GetCl(), 1e-12);
    EXPECT_NEAR(outputs.at("cd")(0), ref.GetCd(), 1e-12);
    EXPECT_NEAR(outputs.at("cm")(0), ref.GetCm(), 1e-12);
}

/* ===========================================================================
   SetOptions error handling
   ===========================================================================
 */

TEST(PhiloteOptionsTest, BothNacaAndCoordinatesThrows) {
    scirocco::SiroccoDiscipline disc;
    disc.Initialize();

    google::protobuf::Struct opts;
    (*opts.mutable_fields())["naca"].set_string_value("0012");
    auto& coord_struct =
        *(*opts.mutable_fields())["coordinates"].mutable_struct_value();
    (*coord_struct.mutable_fields())["x"].mutable_list_value();
    (*coord_struct.mutable_fields())["z"].mutable_list_value();

    EXPECT_THROW(disc.SetOptions(opts), std::invalid_argument);
}

TEST(PhiloteOptionsTest, NeitherNacaNorCoordinatesThrows) {
    scirocco::SiroccoDiscipline disc;
    disc.Initialize();

    google::protobuf::Struct opts;
    (*opts.mutable_fields())["npanel"].set_number_value(99);

    EXPECT_THROW(disc.SetOptions(opts), std::invalid_argument);
}

/* ===========================================================================
   Partial derivative tests: discipline partials vs. CS reference (HDF5)
   ===========================================================================
 */

class PhilotePartialsTest : public ::testing::Test {
   protected:
    static H5::H5File* file_;
    static bool file_loaded_;

    static void SetUpTestSuite() {
        try {
            file_ = new H5::H5File(SENSITIVITY_REFERENCE_FILE, H5F_ACC_RDONLY);
            file_loaded_ = true;
        } catch (H5::Exception&) {
            file_loaded_ = false;
        }
    }

    static void TearDownTestSuite() {
        if (file_) {
            file_->close();
            delete file_;
            file_ = nullptr;
        }
    }

    void SetUp() override {
        if (!file_loaded_)
            GTEST_SKIP() << "Sensitivity reference file not loaded";
    }

    H5::H5File& getFile() { return *file_; }
};

H5::H5File* PhilotePartialsTest::file_ = nullptr;
bool PhilotePartialsTest::file_loaded_ = false;

TEST_F(PhilotePartialsTest, InviscidAlpha) {
    double dcl_ref = readScalar<double>(getFile(), "inviscid_alpha/dcl_dalpha");
    double dcm_ref = readScalar<double>(getFile(), "inviscid_alpha/dcm_dalpha");

    scirocco::SiroccoDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildNacaOptions("2412", 80), inputs, outputs, partials);

    scirocco::Scirocco ref("2412", 80);
    ref.SetOper(5.0);
    ref.Solve();

    SetCoordinateInputs(inputs, ref.GetFoilCoords(), 5.0, 0.0, 0.0);
    disc.Compute(inputs, outputs);
    disc.ComputePartials(inputs, partials);

    auto rel = [](double got, double ref_val) {
        return std::abs(got - ref_val) / std::max(std::abs(ref_val), 1e-10);
    };

    EXPECT_LT(rel(partials[{"cl", "alpha"}](0), dcl_ref), TOL_INVISCID);
    EXPECT_LT(rel(partials[{"cm", "alpha"}](0), dcm_ref), TOL_INVISCID);
}

TEST_F(PhilotePartialsTest, InviscidMach) {
    double dcl_ref = readScalar<double>(getFile(), "inviscid_mach/dcl_dMa");
    double dcm_ref = readScalar<double>(getFile(), "inviscid_mach/dcm_dMa");

    scirocco::SiroccoDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildNacaOptions("0012", 100), inputs, outputs,
                 partials);

    scirocco::Scirocco ref("0012", 100);
    ref.SetOper(3.0, -1.0, 0.0, 0.3);
    ref.Solve();

    SetCoordinateInputs(inputs, ref.GetFoilCoords(), 3.0, 0.0, 0.3);
    disc.Compute(inputs, outputs);
    disc.ComputePartials(inputs, partials);

    auto rel = [](double got, double ref_val) {
        return std::abs(got - ref_val) / std::max(std::abs(ref_val), 1e-10);
    };

    EXPECT_LT(rel(partials[{"cl", "Ma"}](0), dcl_ref), TOL_INVISCID);
    EXPECT_LT(rel(partials[{"cm", "Ma"}](0), dcm_ref), TOL_INVISCID);
}

TEST_F(PhilotePartialsTest, InviscidShape) {
    auto test_nodes = readIntVector(getFile(), "inviscid_shape/test_nodes");
    MatrixXd dcl_dx_ref = readMatrix(getFile(), "inviscid_shape/dcl_dx");
    MatrixXd dcm_dx_ref = readMatrix(getFile(), "inviscid_shape/dcm_dx");
    MatrixXd x_base = readMatrix(getFile(), "inviscid_shape/x_base");

    int N_ref = x_base.cols();
    scirocco::SiroccoDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildNacaOptions("2412", N_ref - 1), inputs, outputs,
                 partials);

    SetCoordinateInputs(inputs, x_base, 5.0, 0.0, 0.0);
    disc.Compute(inputs, outputs);
    disc.ComputePartials(inputs, partials);

    const auto& dcl_dx = partials[{"cl", "x"}];
    const auto& dcm_dx = partials[{"cm", "x"}];

    for (int node : test_nodes) {
        if (node >= N_ref) continue;
        for (int dim = 0; dim < 2; ++dim) {
            int idx = dim * N_ref + node;
            double dcl_got = dcl_dx(idx);
            double dcl_ref_val = dcl_dx_ref(dim, node);
            double rel = std::abs(dcl_got - dcl_ref_val) /
                         std::max(std::abs(dcl_ref_val), 1e-10);
            EXPECT_LT(rel, TOL_INVISCID)
                << "dcl/dx node=" << node << " dim=" << dim;

            double dcm_got = dcm_dx(idx);
            double dcm_ref_val = dcm_dx_ref(dim, node);
            double rel_cm = std::abs(dcm_got - dcm_ref_val) /
                            std::max(std::abs(dcm_ref_val), 1e-10);
            EXPECT_LT(rel_cm, TOL_INVISCID)
                << "dcm/dx node=" << node << " dim=" << dim;
        }
    }
}

TEST_F(PhilotePartialsTest, PythonSplineContract) {
    H5::H5File file(GRADIENT_CONTRACT_REFERENCE_FILE, H5F_ACC_RDONLY);
    for (const std::string name :
         {"inviscid_te", "viscous_spline", "viscous_mach"}) {
        SCOPED_TRACE(name);
        auto group = file.openGroup(name);
        auto attr = [&](const char* key) {
            double value;
            group.openAttribute(key).read(H5::PredType::NATIVE_DOUBLE, &value);
            return value;
        };
        MatrixXd x = readMatrix(file, name + "/x");
        const int n = x.cols();
        scirocco::SiroccoDiscipline disc;
        philote::Variables inputs, outputs;
        philote::Partials partials;
        RunLifecycle(disc, BuildNacaOptions("0012", n - 1), inputs, outputs,
                     partials);
        SetCoordinateInputs(inputs, x, attr("alpha"), attr("Re"), attr("Ma"));
        disc.Compute(inputs, outputs);
        disc.ComputePartials(inputs, partials);
        double primal[3];
        group.openDataSet("outputs").read(primal, H5::PredType::NATIVE_DOUBLE);
        const std::vector<std::string> names = {"cl", "cd", "cm"};
        MatrixXd scalar = readMatrix(file, name + "/scalar_derivatives");
        const std::vector<std::string> params = {"alpha", "Ma", "Re"};
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(outputs.at(names[j])(0), primal[j], 2e-11);
            for (int p = 0; p < 3; ++p) {
                double scale = p == 2 && attr("Re") > 0 ? attr("Re") : 1;
                const double derivative = partials[{names[j], params[p]}](0);
                EXPECT_NEAR(derivative * scale, scalar(p, j) * scale,
                            2e-11 + 2e-8 * std::abs(scalar(p, j) * scale));
            }
        }
        auto data = group.openDataSet("directions");
        hsize_t dims[3];
        data.getSpace().getSimpleExtentDims(dims);
        std::vector<double> directions(dims[0] * 2 * n);
        data.read(directions.data(), H5::PredType::NATIVE_DOUBLE);
        MatrixXd ref = readMatrix(file, name + "/shape_derivatives");
        for (size_t p = 0; p < dims[0]; ++p) {
            SCOPED_TRACE(p);
            for (int j = 0; j < 3; ++j) {
                double got = 0;
                for (int i = 0; i < 2 * n; ++i)
                    got += partials[{names[j], "x"}](i) *
                           directions[p * 2 * n + i];
                EXPECT_NEAR(got, ref(p, j), 2e-11 + 2e-8 * std::abs(ref(p, j)));
            }
        }
        // A previous valid Jacobian must not be reused for changed inputs.
        inputs.at("alpha")(0) += 1;
        EXPECT_THROW(disc.ComputePartials(inputs, partials),
                     std::runtime_error);
    }
}
