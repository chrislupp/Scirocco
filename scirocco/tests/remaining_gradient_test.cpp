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
#include <H5Cpp.h>
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "scirocco/core/adjoint.h"
#include "scirocco/core/paneling.h"
#include "scirocco/core/spline.h"
#include "scirocco/scirocco.h"

using namespace scirocco;
namespace {
MatrixXd ReadMatrix(H5::Group& group, const char* name) {
    auto dataset = group.openDataSet(name);
    hsize_t dims[2];
    dataset.getSpace().getSimpleExtentDims(dims);
    std::vector<double> values(dims[0] * dims[1]);
    dataset.read(values.data(), H5::PredType::NATIVE_DOUBLE);
    MatrixXd out(dims[0], dims[1]);
    for (size_t i = 0; i < dims[0]; ++i)
        for (size_t j = 0; j < dims[1]; ++j)
            out(i, j) = values[i * dims[1] + j];
    return out;
}
double Attribute(H5::Group& group, const char* name) {
    double out;
    group.openAttribute(name).read(H5::PredType::NATIVE_DOUBLE, &out);
    return out;
}
std::vector<MatrixXd> Directions(H5::Group& group) {
    auto data = group.openDataSet("directions");
    hsize_t dims[3];
    data.getSpace().getSimpleExtentDims(dims);
    std::vector<double> values(dims[0] * dims[1] * dims[2]);
    data.read(values.data(), H5::PredType::NATIVE_DOUBLE);
    std::vector<MatrixXd> out;
    for (size_t p = 0; p < dims[0]; ++p) {
        MatrixXd direction(2, dims[2]);
        for (size_t d = 0; d < 2; ++d)
            for (size_t i = 0; i < dims[2]; ++i)
                direction(d, i) = values[(p * 2 + d) * dims[2] + i];
        out.push_back(direction);
    }
    return out;
}
void Check(double got, double reference) {
    EXPECT_TRUE(std::isfinite(got));
    EXPECT_NEAR(got, reference, 2e-11 + 2e-8 * std::abs(reference));
}

class RemainingGradientTest : public ::testing::TestWithParam<const char*> {};
TEST_P(RemainingGradientTest, PythonOutputsAndAllExposedDerivatives) {
    H5::H5File file(REMAINING_GRADIENT_REFERENCE_FILE, H5F_ACC_RDONLY);
    auto group = file.openGroup(GetParam());
    const MatrixXd x = ReadMatrix(group, "x");
    Scirocco foil("0012", x.cols() - 1);
    auto& state = foil.GetState();
    state.foil.x = x;
    if (group.attrExists("chord")) state.geom.chord = Attribute(group, "chord");
    std::string mode;
    auto attr = group.openAttribute("mode");
    attr.read(attr.getStrType(), mode);
    if (mode == "spline")
        RefreshPanelGeometry(state.foil);
    else
        group.openDataSet("s").read(state.foil.s.data(),
                                    H5::PredType::NATIVE_DOUBLE);
    state.oper.redowake = Attribute(group, "redowake") != 0;
    foil.SetOper(Attribute(group, "alpha"), Attribute(group, "target"),
                 Attribute(group, "Re"), Attribute(group, "Ma"));
    ASSERT_TRUE(foil.Solve());
    const Sensitivity sens = foil.ComputeSensitivities();
    if (std::string(GetParam()) == "nonunit_chord") {
        // Independent full Python solve at 80 digits using the original panel
        // closed forms, with central-step and precision convergence checks.
        std::ifstream input(HIGH_PRECISION_GRADIENT_REFERENCE_FILE);
        ASSERT_TRUE(input.good());
        std::string line;
        while (std::getline(input, line) && (line.empty() || line[0] == '#')) {
        }
        std::istringstream fields(line);
        const std::array<double, 4> got{sens.dcl_dx(0, 50), sens.dcd_dx(0, 50),
                                        sens.dcm_dx(0, 50),
                                        sens.dcdpi_dx(0, 50)};
        for (int j = 0; j < 4; ++j) {
            SCOPED_TRACE(j);
            double expected;
            ASSERT_TRUE(static_cast<bool>(fields >> expected));
            EXPECT_NEAR(got[j], expected, 2e-13);
        }
    }
    std::array<double, 5> primal;
    group.openDataSet("outputs").read(primal.data(),
                                      H5::PredType::NATIVE_DOUBLE);
    const std::array<double, 5> outputs{foil.GetCl(), foil.GetCd(),
                                        foil.GetCm(), state.post.cdpi,
                                        state.oper.alpha};
    for (int j = 0; j < 5; ++j) EXPECT_NEAR(outputs[j], primal[j], 3e-10);
    const std::array<MatrixXd, 5> shape{sens.dcl_dx, sens.dcd_dx, sens.dcm_dx,
                                        sens.dcdpi_dx, sens.dalpha_dx};
    const auto directions = Directions(group);
    const MatrixXd reference = ReadMatrix(group, "shape_derivatives");
    for (size_t p = 0; p < directions.size(); ++p) {
        SCOPED_TRACE(p);
        for (int j = 0; j < 5; ++j) {
            SCOPED_TRACE(j);
            Check((shape[j].array() * directions[p].array()).sum(),
                  reference(p, j));
        }
    }
    MatrixXd scalar(4, 5);
    scalar << sens.dcl_dalpha, sens.dcd_dalpha, sens.dcm_dalpha,
        sens.dcdpi_dalpha, sens.dalpha_dalpha, sens.dcl_dMa, sens.dcd_dMa,
        sens.dcm_dMa, sens.dcdpi_dMa, sens.dalpha_dMa, sens.dcl_dRe,
        sens.dcd_dRe, sens.dcm_dRe, sens.dcdpi_dRe, sens.dalpha_dRe,
        sens.dcl_dcltarget, sens.dcd_dcltarget, sens.dcm_dcltarget,
        sens.dcdpi_dcltarget, sens.dalpha_dcltarget;
    const MatrixXd expected = ReadMatrix(group, "scalar_derivatives");
    for (int p = 0; p < 4; ++p)
        for (int j = 0; j < 5; ++j) {
            SCOPED_TRACE(p);
            SCOPED_TRACE(j);
            const double scale =
                p == 2 && state.oper.viscous ? state.oper.Re : 1;
            Check(scalar(p, j) * scale, expected(p, j) * scale);
        }
}
INSTANTIATE_TEST_SUITE_P(
    Mfoil, RemainingGradientTest,
    ::testing::Values("fixed_arc_labels", "negative_alpha", "inviscid_trim",
                      "viscous_trim_frozen", "viscous_trim_rebuilt",
                      "nonunit_chord", "cambered_viscous"));
}  // namespace
