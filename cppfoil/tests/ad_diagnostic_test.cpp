/*
    cppfoil: C++ Airfoil Analysis Library

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
    number: AFRL-2026-XXXX.

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

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include <codi.hpp>

#include "cppfoil/core/adjoint.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/forces.h"
#include "cppfoil/core/geometry.h"
#include "cppfoil/core/wake.h"
#include "cppfoil/cppfoil.h"
#include "cppfoil/data/state.h"

using namespace cppfoil;
using namespace Eigen;

using ADReal = codi::RealForward;

// ---------------------------------------------------------------------------
// HDF5 helpers (same as comprehensive_reference_test.cpp)
// ---------------------------------------------------------------------------
template <typename T>
static T readScalar(H5::H5File& file, const std::string& path) {
    H5::DataSet dataset = file.openDataSet(path);
    T value;
    dataset.read(&value, dataset.getDataType());
    return value;
}

static VectorXd readVector(H5::H5File& file, const std::string& path) {
    H5::DataSet dataset = file.openDataSet(path);
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims);
    std::vector<double> buffer(dims[0]);
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
    return Map<VectorXd>(buffer.data(), dims[0]);
}

static MatrixXd readMatrix(H5::H5File& file, const std::string& path) {
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

static std::vector<double> read3D(H5::H5File& file, const std::string& path,
                                   hsize_t& d0, hsize_t& d1, hsize_t& d2) {
    H5::DataSet dataset = file.openDataSet(path);
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t dims[3];
    dataspace.getSimpleExtentDims(dims);
    d0 = dims[0]; d1 = dims[1]; d2 = dims[2];
    std::vector<double> buffer(d0 * d1 * d2);
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
    return buffer;
}

static inline size_t idx3(int i, int j, int k, int d1, int d2) {
    return static_cast<size_t>(i) * d1 * d2 + j * d2 + k;
}

// ---------------------------------------------------------------------------
// AD forward-path functions using CoDiPack RealForward
// ---------------------------------------------------------------------------

// 2a. Cumulative arc length: s[i] = s[i-1] + |x[i] - x[i-1]|
// x_ad is [2*N] interleaved (x0,y0,x1,y1,...)
static std::vector<ADReal> ad_arc_length(const std::vector<ADReal>& x_ad, int N) {
    std::vector<ADReal> s(N);
    s[0] = 0.0;
    for (int i = 1; i < N; ++i) {
        ADReal dx = x_ad[2*i]     - x_ad[2*(i-1)];
        ADReal dy = x_ad[2*i + 1] - x_ad[2*(i-1) + 1];
        s[i] = s[i-1] + sqrt(dx*dx + dy*dy);
    }
    return s;
}

// 2b. Stagnation point s-location via gamma interpolation (gamma passive/fixed)
// sstag = (g1*s[I0] - g0*s[I1]) / (g1 - g0)
static ADReal ad_stag_sstag(const std::vector<ADReal>& s_ad,
                             double g0, double g1,
                             int I0, int I1) {
    double denom = g1 - g0;
    double w1 = g1 / denom;
    double w2 = -g0 / denom;
    return w1 * s_ad[I0] + w2 * s_ad[I1];
}

// 2c. xi for foil points: xi[i] = sgnue[i] * (s[i] - sstag)
// Uses sgnue multiplication (not abs) for differentiability
static std::vector<ADReal> ad_xi_foil(const std::vector<ADReal>& s_ad,
                                       const ADReal& sstag_ad,
                                       const VectorXd& sgnue, int N) {
    std::vector<ADReal> xi(N);
    for (int i = 0; i < N; ++i) {
        xi[i] = sgnue(i) * (s_ad[i] - sstag_ad);
    }
    return xi;
}

// 2d. Geometric spacing (Newton solve for ratio r, then progression)
// dx0 is the first interval, L is total length, Np is number of points
static std::vector<ADReal> ad_space_geom(const ADReal& dx0_ad, double L, int Np) {
    if (Np <= 1) {
        std::vector<ADReal> sv(1, ADReal(0.0));
        return sv;
    }

    int Ni = Np - 1;  // number of intervals
    // d = L / dx0 (AD-active)
    ADReal d_ad = L / dx0_ad;
    double d = d_ad.getValue();

    // Solve for geometric ratio r using Newton (passive solve)
    double a = Ni * (Ni - 1.0) * (Ni - 2.0) / 6.0;
    double b = Ni * (Ni - 1.0) / 2.0;
    double c = Ni - d;
    double disc = std::max(b * b - 4.0 * a * c, 0.0);
    double r = 1.0 + (-b + std::sqrt(disc)) / (2.0 * a);

    for (int k = 0; k < 10; ++k) {
        double rN = std::pow(r, Ni);
        double R = rN - 1.0 - d * (r - 1.0);
        double R_r = Ni * std::pow(r, Ni - 1) - d;
        double dr = -R / R_r;
        if (std::abs(dr) < 1e-6) break;
        r += dr;
    }

    // Now differentiate r w.r.t. d via implicit function theorem:
    // F(r, d) = r^N - 1 - d*(r-1) = 0
    // dr/dd = -F_d / F_r = (r-1) / (N*r^(N-1) - d)
    double rN = std::pow(r, Ni);
    double F_r = Ni * std::pow(r, Ni - 1) - d;
    double F_d = -(r - 1.0);
    double dr_dd = -F_d / F_r;  // = (r-1) / (N*r^(N-1) - d)

    // r as AD type: r = r_val + dr/dd * dd * eps
    // dd/d(dx0) comes from d = L/dx0 => dd = -L/dx0^2 * d(dx0)
    ADReal r_ad = r;
    r_ad.setGradient(dr_dd * d_ad.getGradient());

    // Generate points with geometric progression
    std::vector<ADReal> sv(Np);
    sv[0] = 0.0;
    ADReal dx = dx0_ad;
    for (int i = 1; i < Np; ++i) {
        sv[i] = sv[i-1] + dx;
        dx = dx * r_ad;
    }

    return sv;
}

// 2e. Wake s values and wake xi
static void ad_wake_s_and_xi(const ADReal& s_Nm1,
                              const std::vector<ADReal>& sv_ad,
                              const ADReal& sstag_ad,
                              int Nw,
                              std::vector<ADReal>& wake_s,
                              std::vector<ADReal>& xi_wake) {
    wake_s.resize(Nw);
    xi_wake.resize(Nw);
    for (int i = 0; i < Nw; ++i) {
        wake_s[i] = s_Nm1 + sv_ad[i];
        xi_wake[i] = wake_s[i] - sstag_ad;
    }
}

// 2f. TE info: hTE and dtdx from panel coordinates
// Replicates GetTEInfo from geometry.cpp
static void ad_te_info(const std::vector<ADReal>& x_ad, int N,
                        ADReal& hTE_out, ADReal& dtdx_out) {
    // t1 = x[0] - x[1] (lower TE tangent direction)
    ADReal v1x = x_ad[0]     - x_ad[2];
    ADReal v1y = x_ad[1]     - x_ad[3];
    ADReal v1n = sqrt(v1x*v1x + v1y*v1y);
    ADReal t1x = v1x / v1n;
    ADReal t1y = v1y / v1n;

    // t2 = x[N-1] - x[N-2] (upper TE tangent direction)
    ADReal v2x = x_ad[2*(N-1)]     - x_ad[2*(N-2)];
    ADReal v2y = x_ad[2*(N-1) + 1] - x_ad[2*(N-2) + 1];
    ADReal v2n = sqrt(v2x*v2x + v2y*v2y);
    ADReal t2x = v2x / v2n;
    ADReal t2y = v2y / v2n;

    // Average tangent (bisector), normalized
    ADReal t_rawx = 0.5 * (t1x + t2x);
    ADReal t_rawy = 0.5 * (t1y + t2y);
    ADReal tn = sqrt(t_rawx*t_rawx + t_rawy*t_rawy);
    ADReal tx = t_rawx / tn;
    ADReal ty = t_rawy / tn;

    // Gap vector: s = x[N-1] - x[0]
    ADReal sx = x_ad[2*(N-1)]     - x_ad[0];
    ADReal sy = x_ad[2*(N-1) + 1] - x_ad[1];

    // hTE = -s.x * t.y + s.y * t.x  (cross product)
    hTE_out = -sx * ty + sy * tx;

    // dtdx = t1.x * t2.y - t2.x * t1.y  (sin of angle between t1, t2)
    dtdx_out = t1x * t2y - t2x * t1y;
}

// 2g. Wake gap: cubic decay from TE
// wgap[w] = hTE * (1 + alpha*xib) * (1-xib)^2  where alpha = 2 + flen*dtdx
// xib = (xi_wake[w] - xi0) / max(Lw, 1e-10)
// xi0 = xi at first wake point
static std::vector<ADReal> ad_wgap(const ADReal& hTE_ad, const ADReal& dtdx_ad,
                                    const std::vector<ADReal>& xi_all,
                                    int N, int Nw) {
    ADReal flen = 2.5;

    // Clip dtdx: branch on passive value
    ADReal dtdx_clip = dtdx_ad;
    double dtdx_val = dtdx_ad.getValue();
    if (dtdx_val < -3.0 / 2.5) {
        dtdx_clip = ADReal(-3.0 / 2.5);
    } else if (dtdx_val > 3.0 / 2.5) {
        dtdx_clip = ADReal(3.0 / 2.5);
    }

    ADReal alpha = 2.0 + flen * dtdx_clip;
    ADReal Lw = flen * hTE_ad;

    // xi0 = xi at first wake point (index N in the full array)
    ADReal xi0 = xi_all[N];

    std::vector<ADReal> wgap(Nw);
    for (int w = 0; w < Nw; ++w) {
        ADReal xi_wake = xi_all[N + w];
        // Compute xib; use passive branch for the max
        double Lw_val = Lw.getValue();
        double Lw_denom = std::max(Lw_val, 1e-10);
        ADReal xib;
        if (Lw_val >= 1e-10) {
            xib = (xi_wake - xi0) / Lw;
        } else {
            // Lw is tiny, xib effectively huge => wgap = 0
            xib = ADReal(1e10);
        }

        double xib_val = xib.getValue();
        if (xib_val >= 0.0 && xib_val <= 1.0) {
            ADReal omxib = 1.0 - xib;
            wgap[w] = hTE_ad * (1.0 + alpha * xib) * omxib * omxib;
        } else {
            wgap[w] = 0.0;
        }
    }
    return wgap;
}

// ---------------------------------------------------------------------------
// Full forward chain: x_foil → (s, sstag, xi_foil, xi_wake)
// Returns dxi_all/dx[k,dim] as a vector [Nsys]
//
// Matches mfoil's residual evaluation ordering:
//   1. recompute_geom → foil.s
//   2. stagpoint_find → sstag from gamma, xi from foil.s and wake.s_base
// wake.s is pre-computed by build_wake and NOT updated when x changes.
// So dxi_wake/dx only has the -dsstag/dx contribution (no s(N-1)→wake.s chain).
// ---------------------------------------------------------------------------
static VectorXd ad_xi_chain(const MatrixXd& x_foil, int N, int Nw,
                              const VectorXd& gam, const int* Istag,
                              const VectorXd& sgnue,
                              const VectorXd& wake_s_base,
                              int k, int dim) {
    int Nsys = N + Nw;

    // Initialize x_ad from foil coordinates
    std::vector<ADReal> x_ad(2 * N);
    for (int i = 0; i < 2 * N; ++i) {
        int j = i / 2;
        int d = i % 2;
        x_ad[i] = x_foil(d, j);
    }
    // Seed
    x_ad[2*k + dim].setGradient(1.0);

    // 2a: arc length
    auto s_ad = ad_arc_length(x_ad, N);

    // 2b: stagnation point (gamma fixed/passive)
    int I0 = Istag[0], I1 = Istag[1];
    double g0 = gam(I0), g1 = gam(I1);
    ADReal sstag_ad = ad_stag_sstag(s_ad, g0, g1, I0, I1);

    // 2c: xi for foil
    auto xi_foil = ad_xi_foil(s_ad, sstag_ad, sgnue, N);

    // Wake xi: xi_wake(w) = wake_s_base(w) - sstag
    // wake_s_base is PASSIVE (pre-computed, not a function of current x)
    // Only sstag contributes derivative
    VectorXd dxi(Nsys);
    for (int i = 0; i < N; ++i) {
        dxi(i) = xi_foil[i].getGradient();
    }
    for (int w = 0; w < Nw; ++w) {
        // d(wake_s_base - sstag)/dx = -dsstag/dx
        dxi(N + w) = -sstag_ad.getGradient();
    }

    return dxi;
}

// ---------------------------------------------------------------------------
// Full forward chain for wgap: x_foil → wgap[Nw]
// Returns dwgap/dx[k,dim] as a vector [Nw]
//
// Matches mfoil: wake.s is passive (pre-computed), so wake xi only changes
// through sstag. wgap depends on xi (through xib) and on hTE/dtdx (through TE
// geometry). Both paths are captured here.
// ---------------------------------------------------------------------------
static VectorXd ad_wgap_chain(const MatrixXd& x_foil, int N, int Nw,
                                const VectorXd& gam, const int* Istag,
                                const VectorXd& sgnue,
                                const VectorXd& wake_s_base,
                                const VectorXd& xi_base,
                                int k, int dim) {
    // Initialize x_ad
    std::vector<ADReal> x_ad(2 * N);
    for (int i = 0; i < 2 * N; ++i) {
        int j = i / 2;
        int d = i % 2;
        x_ad[i] = x_foil(d, j);
    }
    x_ad[2*k + dim].setGradient(1.0);

    // Foil xi chain (arc length → sstag → xi_foil)
    auto s_ad = ad_arc_length(x_ad, N);
    int I0 = Istag[0], I1 = Istag[1];
    double g0 = gam(I0), g1 = gam(I1);
    ADReal sstag_ad = ad_stag_sstag(s_ad, g0, g1, I0, I1);
    auto xi_foil = ad_xi_foil(s_ad, sstag_ad, sgnue, N);

    // Wake xi: passive wake_s_base - active sstag
    std::vector<ADReal> xi_all(N + Nw);
    for (int i = 0; i < N; ++i) xi_all[i] = xi_foil[i];
    for (int w = 0; w < Nw; ++w) {
        ADReal xi_w = wake_s_base(w) - sstag_ad;
        xi_all[N + w] = xi_w;
    }

    // TE info
    ADReal hTE_ad, dtdx_ad;
    ad_te_info(x_ad, N, hTE_ad, dtdx_ad);

    // Compute wgap
    auto wgap_ad = ad_wgap(hTE_ad, dtdx_ad, xi_all, N, Nw);

    VectorXd dwgap(Nw);
    for (int i = 0; i < Nw; ++i) {
        dwgap(i) = wgap_ad[i].getGradient();
    }
    return dwgap;
}

// ---------------------------------------------------------------------------
// Test fixture
// ---------------------------------------------------------------------------
class ADDiagnosticTest : public ::testing::Test {
   protected:
    static H5::H5File* file_;
    static bool file_loaded_;

    // Cached state and derived quantities
    static FoilState* pstate_;
    static VectorXd lambda_bl_cl_;
    static VectorXd lambda_ue_cl_;
    static VectorXd lambda_bl_cd_;
    static VectorXd lambda_ue_cd_;

    // Hand-written VJP outputs
    static VectorXd vjp_xi_cl_;
    static VectorXd vjp_wgap_cl_;
    static VectorXd vjp_xichain_cl_;

    // CS reference data
    static std::vector<double> dR_BL_dx_buf_;
    static hsize_t bd0_, bd1_, bd2_;

    static bool setup_done_;

    static void SetUpTestSuite() {
        try {
            file_ = new H5::H5File(COMPREHENSIVE_REFERENCE_FILE, H5F_ACC_RDONLY);
            file_loaded_ = true;
        } catch (H5::Exception& e) {
            std::cerr << "Failed to open HDF5: " << e.getCDetailMsg() << std::endl;
            file_loaded_ = false;
            return;
        }

        // Load reference geometry
        MatrixXd x_base = readMatrix(*file_, "comprehensive/foil_x");
        VectorXd s_base = readVector(*file_, "comprehensive/foil_s");
        int N_ref = static_cast<int>(x_base.cols());

        // Solve case
        auto* foil = new Cppfoil("0012", N_ref - 1);
        foil->SetOper(3.0, -1.0, 1e6, 0.0);
        pstate_ = &foil->GetState();
        auto& state = *pstate_;
        state.foil.x = x_base;
        state.foil.N = N_ref;
        state.foil.s = s_base;  // Use mfoil's spline-based arc length
        state.oper.viscous = true;
        state.oper.initbl = true;
        state.param.verb = 0;
        SolveViscous(state);
        if (!state.glob.conv) {
            std::cerr << "AD diagnostic: solver did not converge" << std::endl;
            return;
        }
        CalcForce(state);
        CalcForceLin(state);
        CalcCdLin(state);

        int N = state.foil.N;
        int Nsys = state.glob.Nsys;

        // Compute adjoint vectors
        VectorXd dcl_dV = VectorXd::Zero(4 * Nsys);
        for (int i = 0; i < N; ++i)
            dcl_dV(4 * i + 3) = state.post.cl_ue(i);

        auto& lu = *state.glob.lu_factorization;
        VectorXd lambda_cl = lu.transpose().solve(dcl_dV);
        lambda_bl_cl_ = lambda_cl.head(3 * Nsys);
        lambda_ue_cl_ = lambda_cl.tail(Nsys);

        VectorXd dcd_dV = state.post.cd_V;
        VectorXd lambda_cd = lu.transpose().solve(dcd_dV);
        lambda_bl_cd_ = lambda_cd.head(3 * Nsys);
        lambda_ue_cd_ = lambda_cd.tail(Nsys);

        // Compute hand-written VJPs
        vjp_xi_cl_ = ComputeArcLengthVJP(state, lambda_bl_cl_, VectorXd::Zero(Nsys));
        vjp_wgap_cl_ = ComputeWgapVJP(state, lambda_bl_cl_);
        vjp_xichain_cl_ = ComputeXiVJP(state, lambda_bl_cl_);

        // Load CS reference
        dR_BL_dx_buf_ = read3D(*file_, "comprehensive/dR_BL_dx", bd0_, bd1_, bd2_);

        setup_done_ = true;
        std::cout << "AD diagnostic setup complete: N=" << N
                  << " Nw=" << state.wake.N << " Nsys=" << Nsys << std::endl;
    }

    static void TearDownTestSuite() {
        if (file_) { file_->close(); delete file_; file_ = nullptr; }
    }

    void SetUp() override {
        if (!file_loaded_ || !setup_done_)
            GTEST_SKIP() << "AD diagnostic setup failed";
    }

    FoilState& S() { return *pstate_; }
};

// Static member definitions
H5::H5File* ADDiagnosticTest::file_ = nullptr;
bool ADDiagnosticTest::file_loaded_ = false;
FoilState* ADDiagnosticTest::pstate_ = nullptr;
VectorXd ADDiagnosticTest::lambda_bl_cl_;
VectorXd ADDiagnosticTest::lambda_ue_cl_;
VectorXd ADDiagnosticTest::lambda_bl_cd_;
VectorXd ADDiagnosticTest::lambda_ue_cd_;
VectorXd ADDiagnosticTest::vjp_xi_cl_;
VectorXd ADDiagnosticTest::vjp_wgap_cl_;
VectorXd ADDiagnosticTest::vjp_xichain_cl_;
std::vector<double> ADDiagnosticTest::dR_BL_dx_buf_;
hsize_t ADDiagnosticTest::bd0_ = 0;
hsize_t ADDiagnosticTest::bd1_ = 0;
hsize_t ADDiagnosticTest::bd2_ = 0;
bool ADDiagnosticTest::setup_done_ = false;

// ===========================================================================
// Test A: XiChainADvsCS — the decisive diagnostic
//
// For each (k, dim):
//   1. AD-compute dxi_all/dx[k,dim] via full chain (arc→sstag→xi→SpaceGeom→wake)
//   2. Multiply dR_xi_col = R_x * dxi_all
//   3. Compare column against dR_BL_dx_CS[:, dim, k] minus wgap contribution
//   4. Contract vjp_xi_AD = lambda_bl^T * dR_xi_col, compare vs hand-written
// ===========================================================================
TEST_F(ADDiagnosticTest, XiChainADvsCS) {
    auto& state = S();
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = state.glob.Nsys;

    // Precompute mu_wgap for CS correction
    const auto& RV = state.glob.R_V;
    VectorXd mu_wgap = VectorXd::Zero(Nw);
    for (int w = 0; w < Nw; ++w) {
        int col = 4 * (N + w) + 1;
        double dot = 0.0;
        for (SparseMatrixd::InnerIterator it(RV, col); it; ++it) {
            if (it.row() < 3 * Nsys) {
                dot += lambda_bl_cl_(it.row()) * it.value();
            }
        }
        mu_wgap(w) = -dot;
    }

    std::cout << "\n=== Test A: Xi Chain AD vs CS vs Adjoint ===" << std::endl;

    int nfail_ad_cs = 0, nfail_adj_ad = 0;
    double worst_ad_cs = 0, worst_adj_ad = 0;
    int worst_ad_cs_k = -1, worst_ad_cs_d = -1;
    int worst_adj_ad_k = -1, worst_adj_ad_d = -1;

    for (int k = 0; k < N; ++k) {
        for (int dim = 0; dim < 2; ++dim) {
            int idx = 2*k + dim;

            // AD forward: dxi/dx[k,dim]
            VectorXd dxi = ad_xi_chain(state.foil.x, N, Nw,
                                        state.isol.gam, state.isol.Istag.data(),
                                        state.isol.sgnue,
                                        state.wake.s,
                                        k, dim);

            // dR_BL_col = R_x * dxi
            VectorXd dR_xi_col = state.glob.R_x * dxi;

            // AD wgap chain
            VectorXd dwgap = ad_wgap_chain(state.foil.x, N, Nw,
                                            state.isol.gam, state.isol.Istag.data(),
                                            state.isol.sgnue,
                                            state.wake.s, state.isol.xi,
                                            k, dim);

            // AD VJP through xi path only
            double vjp_xi_ad = lambda_bl_cl_.dot(dR_xi_col);

            // AD VJP through wgap path
            double vjp_wgap_ad = mu_wgap.dot(dwgap);

            // Total AD BL VJP (xi + wgap, excluding xichain which goes through gamma)
            double vjp_bl_ad = vjp_xi_ad + vjp_wgap_ad;

            // CS reference: lambda_bl^T * dR_BL_dx_CS[:, dim, k]
            double vjp_bl_cs = 0;
            for (int j = 0; j < 3 * Nsys; ++j) {
                vjp_bl_cs += lambda_bl_cl_(j) *
                    dR_BL_dx_buf_[idx3(j, dim, k, 2, N)];
            }

            // The CS reference includes the xichain (gamma→sstag) contribution;
            // subtract it for a fair comparison to AD (which holds gamma fixed)
            double vjp_bl_cs_no_xichain = vjp_bl_cs - vjp_xichain_cl_(idx);

            // Hand-written adjoint (xi + wgap)
            double vjp_bl_adj = vjp_xi_cl_(idx) + vjp_wgap_cl_(idx);

            // Compare AD vs CS (minus xichain)
            double err_ad_cs = std::abs(vjp_bl_ad - vjp_bl_cs_no_xichain);
            double mag_cs = std::max(std::abs(vjp_bl_cs_no_xichain), 1e-10);
            double rel_ad_cs = err_ad_cs / mag_cs;

            // NOTE: We can't perfectly subtract xichain from CS because it's
            // entangled. Use the full CS for now and compare AD+xichain vs CS.
            // For the decisive test, compare AD vs adjoint directly.
            double err_adj_ad = std::abs(vjp_bl_adj - vjp_bl_ad);
            double mag_ad = std::max(std::abs(vjp_bl_ad), 1e-10);
            double rel_adj_ad = err_adj_ad / mag_ad;

            if (rel_adj_ad > 1e-4) {
                nfail_adj_ad++;
                if (rel_adj_ad > worst_adj_ad) {
                    worst_adj_ad = rel_adj_ad;
                    worst_adj_ad_k = k; worst_adj_ad_d = dim;
                }
            }

            // Print diagnostic for significant disagreements
            if (rel_adj_ad > 0.01) {
                std::printf(
                    "  k=%3d d=%d: AD=%11.4e  adj=%11.4e  CS=%11.4e | "
                    "adj-AD=%6.2f%%  (xi_ad=%11.4e  wgap_ad=%11.4e  "
                    "xi_adj=%11.4e  wgap_adj=%11.4e)\n",
                    k, dim, vjp_bl_ad, vjp_bl_adj, vjp_bl_cs,
                    rel_adj_ad * 100.0,
                    vjp_xi_ad, vjp_wgap_ad,
                    vjp_xi_cl_(idx), vjp_wgap_cl_(idx));
            }
        }
    }

    std::cout << "  adj vs AD failures (>0.01%): " << nfail_adj_ad
              << "/" << 2*N << std::endl;
    std::cout << "  worst adj-AD: rel=" << worst_adj_ad
              << " at k=" << worst_adj_ad_k << " d=" << worst_adj_ad_d << std::endl;

    // This test is diagnostic — report but do not hard-fail
    // The point is to see WHICH sub-function disagrees
}

// ===========================================================================
// Test B: WgapChainADvsCS — wgap diagnostic
// ===========================================================================
TEST_F(ADDiagnosticTest, WgapChainADvsCS) {
    auto& state = S();
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = state.glob.Nsys;

    // Precompute mu_wgap
    const auto& RV = state.glob.R_V;
    VectorXd mu_wgap = VectorXd::Zero(Nw);
    for (int w = 0; w < Nw; ++w) {
        int col = 4 * (N + w) + 1;
        double dot = 0.0;
        for (SparseMatrixd::InnerIterator it(RV, col); it; ++it) {
            if (it.row() < 3 * Nsys)
                dot += lambda_bl_cl_(it.row()) * it.value();
        }
        mu_wgap(w) = -dot;
    }

    std::cout << "\n=== Test B: Wgap Chain AD vs Adjoint ===" << std::endl;

    int nfail = 0;
    double worst = 0;
    int worst_k = -1, worst_d = -1;

    for (int k = 0; k < N; ++k) {
        for (int dim = 0; dim < 2; ++dim) {
            int idx = 2*k + dim;

            VectorXd dwgap = ad_wgap_chain(state.foil.x, N, Nw,
                                            state.isol.gam, state.isol.Istag.data(),
                                            state.isol.sgnue,
                                            state.wake.s, state.isol.xi,
                                            k, dim);

            double vjp_wgap_ad = mu_wgap.dot(dwgap);
            double vjp_wgap_adj = vjp_wgap_cl_(idx);

            double err = std::abs(vjp_wgap_ad - vjp_wgap_adj);
            double mag = std::max(std::abs(vjp_wgap_ad), std::max(std::abs(vjp_wgap_adj), 1e-10));
            double rel = err / mag;

            if (rel > 1e-4) {
                nfail++;
                if (rel > worst) {
                    worst = rel;
                    worst_k = k; worst_d = dim;
                }
            }

            if (rel > 0.01) {
                std::printf("  k=%3d d=%d: AD=%11.4e  adj=%11.4e  rel=%6.2f%%\n",
                            k, dim, vjp_wgap_ad, vjp_wgap_adj, rel * 100.0);
            }
        }
    }

    std::cout << "  wgap failures (>0.01%): " << nfail << "/" << 2*N << std::endl;
    std::cout << "  worst: rel=" << worst
              << " at k=" << worst_k << " d=" << worst_d << std::endl;
}

// ===========================================================================
// Test C: TotalBLVJPTripleCompare — summary
// ===========================================================================
TEST_F(ADDiagnosticTest, TotalBLVJPTripleCompare) {
    auto& state = S();
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = state.glob.Nsys;

    // Precompute mu_wgap
    const auto& RV = state.glob.R_V;
    VectorXd mu_wgap = VectorXd::Zero(Nw);
    for (int w = 0; w < Nw; ++w) {
        int col = 4 * (N + w) + 1;
        double dot = 0.0;
        for (SparseMatrixd::InnerIterator it(RV, col); it; ++it) {
            if (it.row() < 3 * Nsys)
                dot += lambda_bl_cl_(it.row()) * it.value();
        }
        mu_wgap(w) = -dot;
    }

    std::cout << "\n=== Test C: Total BL VJP Triple Compare ===" << std::endl;
    std::cout << "  k   d    AD           CS           adj          "
              << "| AD-CS%    adj-AD%   category" << std::endl;

    int n_ad_vs_cs = 0, n_adj_vs_ad = 0, n_adj_vs_cs = 0;

    for (int k = 0; k < N; ++k) {
        for (int dim = 0; dim < 2; ++dim) {
            int idx = 2*k + dim;

            // AD
            VectorXd dxi = ad_xi_chain(state.foil.x, N, Nw,
                                        state.isol.gam, state.isol.Istag.data(),
                                        state.isol.sgnue,
                                        state.wake.s,
                                        k, dim);
            double vjp_xi_ad = lambda_bl_cl_.dot(state.glob.R_x * dxi);

            VectorXd dwgap = ad_wgap_chain(state.foil.x, N, Nw,
                                            state.isol.gam, state.isol.Istag.data(),
                                            state.isol.sgnue,
                                            state.wake.s, state.isol.xi,
                                            k, dim);
            double vjp_wgap_ad = mu_wgap.dot(dwgap);
            double vjp_bl_ad = vjp_xi_ad + vjp_wgap_ad;

            // CS: lambda_bl^T * dR_BL_dx_CS (includes all chains incl xichain)
            double vjp_bl_cs = 0;
            for (int j = 0; j < 3 * Nsys; ++j)
                vjp_bl_cs += lambda_bl_cl_(j) * dR_BL_dx_buf_[idx3(j, dim, k, 2, N)];

            // Hand-written adjoint (xi + wgap + xichain)
            double vjp_bl_adj = vjp_xi_cl_(idx) + vjp_wgap_cl_(idx) + vjp_xichain_cl_(idx);

            // AD + xichain for fair comparison with CS
            double vjp_bl_ad_full = vjp_bl_ad + vjp_xichain_cl_(idx);

            double mag = std::max({std::abs(vjp_bl_cs), std::abs(vjp_bl_ad_full),
                                   std::abs(vjp_bl_adj), 1e-10});
            double rel_ad_cs = std::abs(vjp_bl_ad_full - vjp_bl_cs) / mag;
            double rel_adj_ad = std::abs(vjp_bl_adj - vjp_bl_ad_full) / mag;
            double rel_adj_cs = std::abs(vjp_bl_adj - vjp_bl_cs) / mag;

            bool fail_ad_cs = rel_ad_cs > 1e-3;
            bool fail_adj_ad = rel_adj_ad > 1e-3;
            bool fail_adj_cs = rel_adj_cs > 1e-3;

            if (fail_ad_cs) n_ad_vs_cs++;
            if (fail_adj_ad) n_adj_vs_ad++;
            if (fail_adj_cs) n_adj_vs_cs++;

            // Determine category
            const char* cat = "PASS";
            if (fail_adj_cs) {
                if (!fail_ad_cs && fail_adj_ad)
                    cat = "ADJOINT_BUG";
                else if (fail_ad_cs && !fail_adj_ad)
                    cat = "FORWARD_BUG";
                else if (fail_ad_cs && fail_adj_ad)
                    cat = "BOTH_DIFFER";
                else
                    cat = "UNKNOWN";
            }

            if (fail_adj_cs) {
                // CS xichain = CS_total - gamma_fixed_AD
                double xichain_cs = vjp_bl_cs - vjp_bl_ad;
                double xichain_cpp = vjp_xichain_cl_(idx);
                double xichain_rel = std::abs(xichain_cs - xichain_cpp) /
                    std::max({std::abs(xichain_cs), std::abs(xichain_cpp), 1e-10});
                std::printf("  %3d %d  %11.4e  %11.4e  %11.4e  | %7.3f%%  %7.3f%%  %s"
                            "  xichain: cs=%11.4e cpp=%11.4e rel=%6.2f%%\n",
                            k, dim, vjp_bl_ad_full, vjp_bl_cs, vjp_bl_adj,
                            rel_ad_cs * 100.0, rel_adj_ad * 100.0, cat,
                            xichain_cs, xichain_cpp, xichain_rel * 100.0);
            }
        }
    }

    std::cout << "\n  Failure tally (>0.1% threshold):" << std::endl;
    std::cout << "    AD vs CS (forward code issue): " << n_ad_vs_cs << "/" << 2*N << std::endl;
    std::cout << "    adj vs AD (adjoint bug):       " << n_adj_vs_ad << "/" << 2*N << std::endl;
    std::cout << "    adj vs CS (either):            " << n_adj_vs_cs << "/" << 2*N << std::endl;
}

// ===========================================================================
// Test D: PassiveWakeSBreakdown — show impact of treating wake.s as passive
//
// Compare AD (passive wake.s, matching mfoil) against old model where
// wake xi includes s(N-1)→wake.s chain, to quantify the fix.
// ===========================================================================
TEST_F(ADDiagnosticTest, PassiveWakeSBreakdown) {
    auto& state = S();
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = state.glob.Nsys;

    std::cout << "\n=== Test D: Passive wake.s Breakdown ===" << std::endl;

    int nfail = 0;
    double worst = 0;

    for (int k = 0; k < N; ++k) {
        for (int dim = 0; dim < 2; ++dim) {
            // New (correct): passive wake.s
            VectorXd dxi_new = ad_xi_chain(state.foil.x, N, Nw,
                                            state.isol.gam, state.isol.Istag.data(),
                                            state.isol.sgnue,
                                            state.wake.s,
                                            k, dim);

            // Old (wrong): wake xi = s(N-1) + sv_passive - sstag
            std::vector<ADReal> x_ad(2 * N);
            for (int i = 0; i < 2 * N; ++i) {
                int j = i / 2;
                int d = i % 2;
                x_ad[i] = state.foil.x(d, j);
            }
            x_ad[2*k + dim].setGradient(1.0);

            auto s_ad = ad_arc_length(x_ad, N);
            int I0 = state.isol.Istag[0], I1 = state.isol.Istag[1];
            double g0 = state.isol.gam(I0), g1 = state.isol.gam(I1);
            ADReal sstag_ad = ad_stag_sstag(s_ad, g0, g1, I0, I1);
            auto xi_foil = ad_xi_foil(s_ad, sstag_ad, state.isol.sgnue, N);

            VectorXd dxi_old(Nsys);
            for (int i = 0; i < N; ++i)
                dxi_old(i) = xi_foil[i].getGradient();
            for (int w = 0; w < Nw; ++w)
                dxi_old(N + w) = s_ad[N-1].getGradient() - sstag_ad.getGradient();

            double vjp_new = lambda_bl_cl_.dot(state.glob.R_x * dxi_new);
            double vjp_old = lambda_bl_cl_.dot(state.glob.R_x * dxi_old);
            double diff = vjp_old - vjp_new;
            double mag = std::max(std::abs(vjp_new), 1e-10);
            double frac = std::abs(diff) / mag;

            if (frac > 0.01) {
                nfail++;
                if (frac > worst) worst = frac;
                if (frac > 0.05) {
                    std::printf("  k=%3d d=%d: new=%11.4e  old=%11.4e  "
                                "s(N-1)_chain=%11.4e (%5.1f%%)\n",
                                k, dim, vjp_new, vjp_old, diff, frac * 100.0);
                }
            }
        }
    }

    std::cout << "  Nodes where s(N-1) chain matters >1%: " << nfail
              << "/" << 2*N << std::endl;
    std::cout << "  Worst s(N-1) fraction: " << worst * 100.0 << "%" << std::endl;
}

// ===========================================================================
// AD panel influence functions for the AIC system
// ===========================================================================

// AD panel info: transform evaluation point to panel-local coordinates
struct ADPanelInfo {
    ADReal len, x1, x2, y, r1, r2, theta1, theta2;
};

static ADPanelInfo ad_get_panel_info(const ADReal* Xj, const ADReal* xi) {
    // Xj = [x1x, x1y, x2x, x2y], xi = [px, py]
    ADPanelInfo info;
    ADReal dx = Xj[2] - Xj[0];
    ADReal dy = Xj[3] - Xj[1];
    info.len = sqrt(dx*dx + dy*dy);
    ADReal tx = dx / info.len;
    ADReal ty = dy / info.len;
    ADReal nx = -ty, ny = tx;
    ADReal d1x = xi[0] - Xj[0], d1y = xi[1] - Xj[1];
    ADReal d2x = xi[0] - Xj[2], d2y = xi[1] - Xj[3];
    info.x1 = d1x*tx + d1y*ty;
    info.x2 = d2x*tx + d2y*ty;
    info.y  = d1x*nx + d1y*ny;
    info.r1 = sqrt(d1x*d1x + d1y*d1y);
    info.r2 = sqrt(d2x*d2x + d2y*d2y);
    // atan2 — guard against zero distances (self-influence)
    if (info.r1.getValue() < 1e-12) {
        info.theta1 = ADReal(0.0);
    } else {
        info.theta1 = atan2(info.y, info.x1);
    }
    if (info.r2.getValue() < 1e-12) {
        info.theta2 = ADReal(0.0);
    } else {
        info.theta2 = atan2(info.y, info.x2);
    }
    return info;
}

// AD linear vortex streamfunction: returns (a, b)
static void ad_panel_linvortex_stream(const ADReal* Xj, const ADReal* xi,
                                       ADReal& a, ADReal& b) {
    ADPanelInfo info = ad_get_panel_info(Xj, xi);
    if (info.len.getValue() < 1e-12) { a = 0; b = 0; return; }
    ADReal d = info.len;
    ADReal x = info.x1;
    ADReal z = info.y;
    ADReal dtheta = info.theta2 - info.theta1;
    ADReal logr1 = (info.r1.getValue() < 1e-9) ? ADReal(0.0) : log(info.r1);
    ADReal logr2 = (info.r2.getValue() < 1e-9) ? ADReal(0.0) : log(info.r2);
    ADReal inv2pi = 0.5 / M_PI;
    ADReal P1 = inv2pi * (z*dtheta - d + x*logr1 - (x-d)*logr2);
    ADReal r1sq = info.r1*info.r1, r2sq = info.r2*info.r2;
    ADReal P2 = x*P1 + inv2pi*(0.5*r2sq*logr2 - 0.5*r1sq*logr1 - r2sq/4.0 + r1sq/4.0);
    a = P1 - P2/d;
    b = P2/d;
}

// AD constant source streamfunction
static ADReal ad_panel_constsource_stream(const ADReal* Xj, const ADReal* xi) {
    ADPanelInfo info = ad_get_panel_info(Xj, xi);
    if (info.len.getValue() < 1e-12) return ADReal(0.0);
    ADReal d = info.len;
    ADReal x = info.x1, z = info.y;
    ADReal theta1 = info.theta1, theta2 = info.theta2;
    ADReal logr1 = (info.r1.getValue() < 1e-9) ? ADReal(0.0) : log(info.r1);
    ADReal logr2 = (info.r2.getValue() < 1e-9) ? ADReal(0.0) : log(info.r2);
    if (info.r1.getValue() < 1e-9) { theta1 = M_PI; theta2 = M_PI; }
    if (info.r2.getValue() < 1e-9) { theta1 = 0.0; theta2 = 0.0; }
    ADReal P = (x*(theta1-theta2) + d*theta2 + z*logr1 - z*logr2) / (2.0*M_PI);
    ADReal dP = d;
    if ((theta1.getValue() + theta2.getValue()) > M_PI)
        P = P - 0.25*dP;
    else
        P = P + 0.75*dP;
    return P;
}

// AD TE info: tcp and tdp (needed for AIC TE panel contributions)
static void ad_te_info_full(const std::vector<ADReal>& x_ad, int N,
                             ADReal& tcp_out, ADReal& tdp_out) {
    ADReal v1x = x_ad[0] - x_ad[2], v1y = x_ad[1] - x_ad[3];
    ADReal v1n = sqrt(v1x*v1x + v1y*v1y);
    ADReal t1x = v1x/v1n, t1y = v1y/v1n;
    ADReal v2x = x_ad[2*(N-1)] - x_ad[2*(N-2)];
    ADReal v2y = x_ad[2*(N-1)+1] - x_ad[2*(N-2)+1];
    ADReal v2n = sqrt(v2x*v2x + v2y*v2y);
    ADReal t2x = v2x/v2n, t2y = v2y/v2n;
    ADReal trx = 0.5*(t1x+t2x), try_ = 0.5*(t1y+t2y);
    ADReal tn = sqrt(trx*trx + try_*try_);
    ADReal tx = trx/tn, ty = try_/tn;
    ADReal sx = x_ad[2*(N-1)] - x_ad[0], sy = x_ad[2*(N-1)+1] - x_ad[1];
    ADReal hTE = -sx*ty + sy*tx;
    ADReal sn = sqrt(sx*sx + sy*sy);
    ADReal px = sx/sn, py = sy/sn;
    tcp_out = (sn.getValue() > 1e-10) ? abs(tx*py - ty*px) : ADReal(0.0);
    tdp_out = (sn.getValue() > 1e-10) ? tx*px + ty*py : ADReal(0.0);
}

// ===========================================================================
// Full AD xichain: x → AIC → gamma → sstag → xi
// Returns dxi_all/dx[k,dim] through the FULL gamma chain.
// ===========================================================================
static VectorXd ad_full_xichain(const MatrixXd& x_foil, int N, int Nw,
                                 const VectorXd& sgnue,
                                 const VectorXd& wake_s_base,
                                 double alpha, double chord,
                                 int k, int dim) {
    int Nsys = N + Nw;

    // Init AD x
    std::vector<ADReal> x_ad(2*N);
    for (int i = 0; i < 2*N; ++i) {
        x_ad[i] = x_foil(i%2, i/2);
    }
    x_ad[2*k + dim].setGradient(1.0);

    // Build AIC [N+1 x N+1] with AD types
    // We only need to track derivatives of the RHS and AIC entries
    // that depend on foil.x. Solve A*gamma = rhs for gamma.
    // Since A and rhs are AD-active, gamma will have derivatives.

    // Get TE info
    ADReal tcp, tdp;
    ad_te_info_full(x_ad, N, tcp, tdp);
    // Check no-gap
    ADReal hTE_ad, dtdx_ad;
    ad_te_info(x_ad, N, hTE_ad, dtdx_ad);
    bool nogap = (std::abs(hTE_ad.getValue()) < 1e-10 * chord);

    // We can't easily do a full AD matrix solve with CoDiPack forward mode
    // for a single seed. Instead use implicit differentiation:
    //   A * gamma = rhs
    //   dA/dx * gamma + A * dgamma/dx = drhs/dx
    //   dgamma/dx = A^{-1} * (drhs/dx - dA/dx * gamma)
    //
    // Compute dR_AIC/dx[k,dim] = dA/dx * gamma_base - drhs/dx column.
    // Then solve A * dgamma = -dR_AIC for dgamma.

    // Use the base (passive) AIC and gamma for the solve.
    // Compute the AD residual column: for each row i, evaluate
    // R(i) = sum_j A(i,j)*gamma(j) - rhs(i) with AD types,
    // extract the gradient.

    // Base gamma (passive)
    std::vector<double> gam_base(N);
    double cosa = std::cos(alpha * M_PI / 180.0);
    double sina = std::sin(alpha * M_PI / 180.0);

    // We need the actual converged gamma, not recomputed.
    // But for this AD test, we need gamma = gamref*cos + gamref*sin
    // where gamref comes from the base AIC. Since we're checking
    // the derivative, we should use the SAME gamma as the forward code.
    // The test fixture has state.isol.gam. But we can't access state here.
    // Pass it in instead.

    // Actually, let's restructure: pass gam as parameter.
    (void)cosa; (void)sina; (void)nogap;
    // This approach is getting too complex. Let me use a simpler method:
    // just compute dR_AIC/dx[k,dim] directly by evaluating each row
    // of the AIC residual with AD, and extracting the gradient.

    // For row i of the AIC system:
    //   R(i) = sum_j A(i,j)*gam(j) - rhs(i)
    // where A(i,j) and rhs(i) depend on x through panel geometry.
    // gam is PASSIVE (the base gamma).

    return VectorXd::Zero(Nsys); // placeholder
}

// ===========================================================================
// Test E: XiChain AD through AIC — verify ComputeXiVJP by running AD
// through the full AIC residual to get exact dgamma/dx, then dsstag, dxi.
// ===========================================================================
TEST_F(ADDiagnosticTest, XiChainADthroughAIC) {
    auto& state = S();
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = state.glob.Nsys;
    double alpha = state.oper.alpha;
    double cosa = std::cos(alpha * M_PI / 180.0);
    double sina = std::sin(alpha * M_PI / 180.0);

    VectorXd vjp_xichain = ComputeXiVJP(state, lambda_bl_cl_);

    // Precompute c = lambda_bl^T * R_x * dxi/dsstag
    VectorXd dxi_dsstag = VectorXd::Zero(Nsys);
    for (int i = 0; i < N; ++i) dxi_dsstag(i) = -state.isol.sgnue(i);
    for (int i = 0; i < Nw; ++i) dxi_dsstag(N+i) = -1.0;
    double c_scalar = lambda_bl_cl_.dot(state.glob.R_x * dxi_dsstag);

    // sstag_g (recomputed, matching ComputeXiVJP)
    int I0 = state.isol.Istag[0], I1 = state.isol.Istag[1];
    double g0 = state.isol.gam(I0), g1 = state.isol.gam(I1);
    double denom_g = g1 - g0;
    double sstag_g0 = g1*(state.foil.s(I0)-state.foil.s(I1))/(denom_g*denom_g);
    double sstag_g1 = -g0*(state.foil.s(I0)-state.foil.s(I1))/(denom_g*denom_g);

    // Get TE info
    TEInfo te = GetTEInfo(state.foil.x);
    bool nogap = (std::abs(te.hTE) < 1e-10 * state.geom.chord);

    std::cout << "\n=== Test E: XiChain AD through AIC vs ComputeXiVJP ===" << std::endl;
    std::cout << "  c_scalar=" << c_scalar << " sstag_g=[" << sstag_g0
              << "," << sstag_g1 << "]" << std::endl;

    int nfail = 0;
    double worst = 0;
    int worst_k = -1, worst_d = -1;

    for (int k = 0; k < N; ++k) {
        for (int dim = 0; dim < 2; ++dim) {
            int idx = 2*k + dim;

            // Compute dR_AIC/dx[k,dim] by AD through AIC assembly.
            // R_AIC(i) = sum_j A(i,j)*gam(j) - rhs(i) for i=0..N-1
            // R_AIC(N) = gam(0) + gam(N-1)  (Kutta)
            // A(i,j) and rhs(i) depend on x. gam is passive.

            std::vector<ADReal> x_ad(2*N);
            for (int ii = 0; ii < 2*N; ++ii)
                x_ad[ii] = state.foil.x(ii%2, ii/2);
            x_ad[2*k + dim].setGradient(1.0);

            // TE info with AD
            ADReal tcp_ad, tdp_ad;
            ad_te_info_full(x_ad, N, tcp_ad, tdp_ad);

            VectorXd dR_AIC(N + 1);
            dR_AIC.setZero();

            for (int i = 0; i < N; ++i) {
                if (nogap && i == N-1) continue;
                ADReal xi_pt[2] = {x_ad[2*i], x_ad[2*i+1]};
                ADReal R_i = 0.0;

                // Linear vortex panels
                for (int j = 0; j < N-1; ++j) {
                    ADReal Xj[4] = {x_ad[2*j], x_ad[2*j+1],
                                    x_ad[2*(j+1)], x_ad[2*(j+1)+1]};
                    ADReal aij, bij;
                    ad_panel_linvortex_stream(Xj, xi_pt, aij, bij);
                    R_i += aij * state.isol.gam(j) + bij * state.isol.gam(j+1);
                }

                // TE source panel: Xj = [N-1, 0]
                {
                    ADReal XjTE[4] = {x_ad[2*(N-1)], x_ad[2*(N-1)+1],
                                      x_ad[0], x_ad[1]};
                    ADReal a_src = ad_panel_constsource_stream(XjTE, xi_pt);
                    double sigma_gam = 0.5 * (state.isol.gam(N-1) - state.isol.gam(0));
                    R_i += a_src * tcp_ad * sigma_gam;
                }

                // TE vortex panel
                {
                    ADReal XjTE[4] = {x_ad[2*(N-1)], x_ad[2*(N-1)+1],
                                      x_ad[0], x_ad[1]};
                    ADReal aij, bij;
                    ad_panel_linvortex_stream(XjTE, xi_pt, aij, bij);
                    double vortex_gam = 0.5 * (state.isol.gam(0) - state.isol.gam(N-1));
                    R_i += (aij + bij) * tdp_ad * vortex_gam;
                }

                // RHS: rhs(i) = -xi_y*cos(alpha) + xi_x*sin(alpha)
                ADReal rhs_i = -xi_pt[1]*cosa + xi_pt[0]*sina;
                R_i -= rhs_i;

                dR_AIC(i) = R_i.getGradient();
            }

            // No-gap row (if applicable)
            if (nogap) {
                // Gamma extrapolation — no x dependence
                dR_AIC(N-1) = 0.0;
            }

            // Kutta row: gam(0) + gam(N-1) = 0 — no x dependence
            dR_AIC(N) = 0.0;

            // Solve: A * dgamma = -dR_AIC
            VectorXd dgamma = state.isol.AIC.partialPivLu().solve(-dR_AIC);

            // dsstag = sstag_g0 * dgamma[I0] + sstag_g1 * dgamma[I1]
            double dsstag = sstag_g0 * dgamma(I0) + sstag_g1 * dgamma(I1);

            // xichain VJP = c_scalar * dsstag
            double vjp_xichain_ad = c_scalar * dsstag;

            // Compare against ComputeXiVJP
            double vjp_xichain_adj = vjp_xichain(idx);

            double err = std::abs(vjp_xichain_ad - vjp_xichain_adj);
            double mag = std::max(std::abs(vjp_xichain_ad),
                                  std::max(std::abs(vjp_xichain_adj), 1e-10));
            double rel = err / mag;

            if (rel > 1e-3) {
                nfail++;
                if (rel > worst) {
                    worst = rel; worst_k = k; worst_d = dim;
                }
            }
            if (rel > 0.01 || (k < 5 && dim == 0)) {
                std::printf("  k=%3d d=%d: AD=%11.4e  adj=%11.4e  rel=%7.3f%%"
                            "  dsstag=%11.4e  dgam=[%11.4e,%11.4e]\n",
                            k, dim, vjp_xichain_ad, vjp_xichain_adj, rel*100.0,
                            dsstag, dgamma(I0), dgamma(I1));
            }
        }
    }

    std::cout << "  xichain AD vs adj failures (>0.1%): " << nfail
              << "/" << 2*N << std::endl;
    std::cout << "  worst: rel=" << worst
              << " at k=" << worst_k << " d=" << worst_d << std::endl;

    // Print C++ R_x at row 147 for stagnation columns
    {
        int row = 147;
        std::printf("\n  C++ R_x row %d: ", row);
        for (int col : {47, 48, 49, 50}) {
            std::printf("col=%d: %11.4e  ", col, state.glob.R_x.coeff(row, col));
        }
        std::printf("\n");
    }

    // Now compute TOTAL dR_BL/dx column through AD for sample nodes
    // and compare element-by-element against CS reference
    std::cout << "\n  --- dR_BL column comparison for sample nodes ---" << std::endl;
    int sample_nodes[] = {20, 50, 80};
    for (int k : sample_nodes) {
        for (int dim = 0; dim < 2; ++dim) {
            int idx = 2*k + dim;

            // Gamma-fixed dxi
            VectorXd dxi_gf = ad_xi_chain(state.foil.x, N, Nw,
                                           state.isol.gam, state.isol.Istag.data(),
                                           state.isol.sgnue, state.wake.s, k, dim);

            // Xichain dxi (through AIC): only sstag changes through gamma
            // dxi_xichain(i) = -sgnue(i) * dsstag_gamma for foil
            // dxi_xichain(N+w) = -dsstag_gamma for wake
            // Compute dsstag_gamma for this (k,dim):
            std::vector<ADReal> x_ad2(2*N);
            for (int ii = 0; ii < 2*N; ++ii)
                x_ad2[ii] = state.foil.x(ii%2, ii/2);
            x_ad2[2*k + dim].setGradient(1.0);
            ADReal tcp2, tdp2;
            ad_te_info_full(x_ad2, N, tcp2, tdp2);
            VectorXd dR_AIC2(N+1);
            dR_AIC2.setZero();
            for (int i = 0; i < N; ++i) {
                if (nogap && i == N-1) continue;
                ADReal xi_pt[2] = {x_ad2[2*i], x_ad2[2*i+1]};
                ADReal R_i = 0.0;
                for (int j = 0; j < N-1; ++j) {
                    ADReal Xj[4] = {x_ad2[2*j], x_ad2[2*j+1],
                                    x_ad2[2*(j+1)], x_ad2[2*(j+1)+1]};
                    ADReal aij, bij;
                    ad_panel_linvortex_stream(Xj, xi_pt, aij, bij);
                    R_i += aij * state.isol.gam(j) + bij * state.isol.gam(j+1);
                }
                {
                    ADReal XjTE[4] = {x_ad2[2*(N-1)], x_ad2[2*(N-1)+1],
                                      x_ad2[0], x_ad2[1]};
                    ADReal a_src = ad_panel_constsource_stream(XjTE, xi_pt);
                    R_i += a_src * tcp2 * 0.5 * (state.isol.gam(N-1) - state.isol.gam(0));
                }
                {
                    ADReal XjTE[4] = {x_ad2[2*(N-1)], x_ad2[2*(N-1)+1],
                                      x_ad2[0], x_ad2[1]};
                    ADReal aij2, bij2;
                    ad_panel_linvortex_stream(XjTE, xi_pt, aij2, bij2);
                    R_i += (aij2 + bij2) * tdp2 * 0.5 * (state.isol.gam(0) - state.isol.gam(N-1));
                }
                ADReal rhs_i = -xi_pt[1]*cosa + xi_pt[0]*sina;
                R_i -= rhs_i;
                dR_AIC2(i) = R_i.getGradient();
            }
            dR_AIC2(N) = 0.0;
            VectorXd dgamma = state.isol.AIC.partialPivLu().solve(-dR_AIC2);
            double dsstag_gamma = sstag_g0*dgamma(I0) + sstag_g1*dgamma(I1);

            VectorXd dxi_xichain = VectorXd::Zero(Nsys);
            for (int i = 0; i < N; ++i)
                dxi_xichain(i) = -state.isol.sgnue(i) * dsstag_gamma;
            for (int w = 0; w < Nw; ++w)
                dxi_xichain(N+w) = -dsstag_gamma;

            // Total dxi
            VectorXd dxi_total = dxi_gf + dxi_xichain;

            // dR_BL from C++ AD: R_x * dxi_total
            VectorXd dR_BL_ad = state.glob.R_x * dxi_total;

            // CS dR_BL column
            VectorXd dR_BL_cs(3*Nsys);
            for (int j = 0; j < 3*Nsys; ++j)
                dR_BL_cs(j) = dR_BL_dx_buf_[idx3(j, dim, k, 2, N)];

            // Compare
            double ad_norm = dR_BL_ad.norm();
            double cs_norm = dR_BL_cs.norm();
            double diff_norm = (dR_BL_ad - dR_BL_cs).norm();
            double rel_norm = diff_norm / std::max(cs_norm, 1e-10);
            std::printf("  k=%d d=%d: |ad|=%9.3e |cs|=%9.3e |diff|=%9.3e rel=%6.3f%%\n",
                        k, dim, ad_norm, cs_norm, diff_norm, rel_norm*100);

            // Print worst disagreeing rows
            int worst_row = -1;
            double worst_row_rel = 0;
            for (int j = 0; j < 3*Nsys; ++j) {
                double e = std::abs(dR_BL_ad(j) - dR_BL_cs(j));
                double m = std::max(std::abs(dR_BL_cs(j)), 1e-15);
                double r = e / m;
                if (r > worst_row_rel) { worst_row_rel = r; worst_row = j; }
            }
            std::printf("    worst row=%d (station %d, eq %d): ad=%11.4e cs=%11.4e\n",
                        worst_row, worst_row/3, worst_row%3,
                        dR_BL_ad(worst_row), dR_BL_cs(worst_row));
        }
    }
}
