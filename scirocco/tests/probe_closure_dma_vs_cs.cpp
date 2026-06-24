/**
 * @brief Compare scirocco's analytic _dMa helpers against mfoil_cs complex-step
 * reference (closure_dma_cs_reference.h5). No FD anywhere.
 */

#include <H5Cpp.h>

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

#include "scirocco/core/bl_chain_internals.h"
#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/properties/bl_properties.h"
#include "scirocco/scirocco.h"

using namespace scirocco;
using namespace scirocco::detail;
using namespace Eigen;

static VectorXd readVec(H5::H5File& f, const std::string& p) {
    H5::DataSet d = f.openDataSet(p);
    hsize_t dims[1];
    d.getSpace().getSimpleExtentDims(dims);
    std::vector<double> b(dims[0]);
    d.read(b.data(), H5::PredType::NATIVE_DOUBLE);
    return Map<VectorXd>(b.data(), dims[0]);
}
static std::vector<int8_t> readI8(H5::H5File& f, const std::string& p) {
    H5::DataSet d = f.openDataSet(p);
    hsize_t dims[1];
    d.getSpace().getSimpleExtentDims(dims);
    std::vector<int8_t> b(dims[0]);
    d.read(b.data(), H5::PredType::NATIVE_INT8);
    return b;
}
static MatrixXd readMat(H5::H5File& f, const std::string& p) {
    H5::DataSet d = f.openDataSet(p);
    hsize_t dims[2];
    d.getSpace().getSimpleExtentDims(dims);
    std::vector<double> b(dims[0] * dims[1]);
    d.read(b.data(), H5::PredType::NATIVE_DOUBLE);
    MatrixXd m(dims[0], dims[1]);
    for (hsize_t i = 0; i < dims[0]; ++i)
        for (hsize_t j = 0; j < dims[1]; ++j) m(i, j) = b[i * dims[1] + j];
    return m;
}

/**
 * @brief Entry point: loads reference data, runs analytic _dMa helpers, and
 *        reports max absolute and relative errors per closure.
 * @return 0 on success.
 */
int main() {
    const std::string path =
        std::string(SCIROCCO_TEST_DATA_DIR) + "/closure_dma_cs_reference.h5";
    H5::H5File f(path, H5F_ACC_RDONLY);
    int N = (int)readVec(f, "N")(0);
    int Nw = (int)readVec(f, "Nw")(0);
    int Nsys = N + Nw;
    double Ma_ref = readVec(f, "Ma")(0);
    MatrixXd U_ref = readMat(f, "U");
    VectorXd xi_ref = readVec(f, "xi");
    auto turb_flags = readI8(f, "turb");
    VectorXd cs_GetUk = readVec(f, "GetUk");
    VectorXd cs_GetMach2 = readVec(f, "GetMach2");
    VectorXd cs_GetHk = readVec(f, "GetHk");
    VectorXd cs_GetHs_T = readVec(f, "GetHs_T");
    VectorXd cs_GetHs_F = readVec(f, "GetHs_F");
    VectorXd cs_GetRet = readVec(f, "GetRet");
    VectorXd cs_GetUs_T = readVec(f, "GetUs_T");
    VectorXd cs_GetUs_F = readVec(f, "GetUs_F");
    VectorXd cs_GetCf_T = readVec(f, "GetCf_T");
    VectorXd cs_GetCf_F = readVec(f, "GetCf_F");
    VectorXd cs_GetCDi_T = readVec(f, "GetCDi_T");
    VectorXd cs_GetCDi_F = readVec(f, "GetCDi_F");
    VectorXd cs_GetCfxt_T = readVec(f, "GetCfxt_T");
    VectorXd cs_GetCfxt_F = readVec(f, "GetCfxt_F");
    VectorXd cs_GetCDixt_T = readVec(f, "GetCDixt_T");
    VectorXd cs_GetCDixt_F = readVec(f, "GetCDixt_F");

    Scirocco foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, Ma_ref);
    auto& s = foil.GetState();
    s.oper.viscous = true;
    s.oper.initbl = true;
    s.param.verb = 0;
    SolveViscous(s);

    /* Use the scirocco-side U values (they should match mfoil's U_ref at
     * converged state for the same setup). */
    Param base_p = s.param;
    ThermoMaDeriv d = ComputeThermoMaDeriv(base_p, base_p.Vinf, 0.0);

    auto report = [](const std::string& name, const VectorXd& cs,
                     const std::vector<double>& an) {
        double max_abs = 0, max_rel = 0;
        int wi = 0;
        for (int i = 0; i < (int)an.size(); ++i) {
            double diff = std::abs(cs(i) - an[i]);
            double mag = std::max(std::abs(cs(i)), std::abs(an[i]));
            if (diff > max_abs) {
                max_abs = diff;
                wi = i;
            }
            if (mag > 1e-8) {
                double rel = diff / mag;
                if (rel > max_rel) max_rel = rel;
            }
        }
        std::cout << "  " << name << ": max abs=" << max_abs << " (at i=" << wi
                  << " cs=" << cs(wi) << " an=" << an[wi]
                  << ") max rel=" << max_rel << "\n";
    };

    std::cout << "Closure ∂/∂Ma: scirocco analytic vs mfoil CS reference\n\n";

    int Nsys_local = s.glob.Nsys;
    std::vector<double> an_Uk(Nsys), an_M2(Nsys), an_Hk(Nsys);
    std::vector<double> an_Hs_T(Nsys), an_Hs_F(Nsys), an_Ret(Nsys);
    std::vector<double> an_Us_T(Nsys), an_Us_F(Nsys);
    std::vector<double> an_Cf_T(Nsys), an_Cf_F(Nsys);
    std::vector<double> an_CDi_T(Nsys), an_CDi_F(Nsys);
    std::vector<double> an_Cfxt_T(Nsys), an_Cfxt_F(Nsys);
    std::vector<double> an_CDixt_T(Nsys), an_CDixt_F(Nsys);

    for (int node = 0; node < Nsys; ++node) {
        Vector4d U = s.glob.U.col(node);
        bool is_wake = (node >= N);
        bool is_turb = (node < N) ? (bool)turb_flags[node] : true;
        double xi = s.isol.xi(node);

        auto make_p = [&](bool turb_flag, bool wake_flag) {
            Param p = base_p;
            p.turb = turb_flag;
            p.wake = wake_flag;
            p.simi = false;
            return p;
        };
        Param p_nat = make_p(is_turb, is_wake);
        Param p_t = make_p(true, is_wake);
        Param p_f = make_p(false, is_wake);

        an_Uk[node] = GetUk_dMa(U(3), p_nat, d);
        an_M2[node] = GetMach2_dMa(U, p_nat, d);
        an_Hk[node] = GetHk_dMa(U, p_nat, d);
        an_Hs_T[node] = GetHs_dMa(U, p_t, d);
        an_Hs_F[node] = GetHs_dMa(U, p_f, d);
        an_Ret[node] = GetRet_dMa(U, p_nat, d);
        an_Us_T[node] = GetUs_dMa(U, p_t, d);
        an_Us_F[node] = GetUs_dMa(U, p_f, d);
        an_Cf_T[node] = GetCf_dMa(U, p_t, d);
        an_Cf_F[node] = GetCf_dMa(U, p_f, d);
        an_CDi_T[node] = GetCDi_dMa(U, p_t, d);
        an_CDi_F[node] = GetCDi_dMa(U, p_f, d);
        an_Cfxt_T[node] = GetCfxt_dMa(U, xi, p_t, d);
        an_Cfxt_F[node] = GetCfxt_dMa(U, xi, p_f, d);
        an_CDixt_T[node] = GetCDixt_dMa(U, xi, p_t, d);
        an_CDixt_F[node] = GetCDixt_dMa(U, xi, p_f, d);
    }

    report("GetUk     ", cs_GetUk, an_Uk);
    report("GetMach2  ", cs_GetMach2, an_M2);
    report("GetHk     ", cs_GetHk, an_Hk);
    report("GetHs(T)  ", cs_GetHs_T, an_Hs_T);
    report("GetHs(F)  ", cs_GetHs_F, an_Hs_F);
    report("GetRet    ", cs_GetRet, an_Ret);
    report("GetUs(T)  ", cs_GetUs_T, an_Us_T);
    report("GetUs(F)  ", cs_GetUs_F, an_Us_F);
    report("GetCf(T)  ", cs_GetCf_T, an_Cf_T);
    report("GetCf(F)  ", cs_GetCf_F, an_Cf_F);
    report("GetCDi(T) ", cs_GetCDi_T, an_CDi_T);
    report("GetCDi(F) ", cs_GetCDi_F, an_CDi_F);
    report("GetCfxt(T)", cs_GetCfxt_T, an_Cfxt_T);
    report("GetCfxt(F)", cs_GetCfxt_F, an_Cfxt_F);
    report("GetCDixt(T)", cs_GetCDixt_T, an_CDixt_T);
    report("GetCDixt(F)", cs_GetCDixt_F, an_CDixt_F);

    (void)Nsys_local;
    return 0;
}
