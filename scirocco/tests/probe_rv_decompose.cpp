/**
 * @brief Decompose R_V[144, 195] (R_mom@node48 wrt U[ue, node48]) into:
 *   (a) the simi-station R1_U(0, 3) contribution from StagnationState chain
 *   (b) the xi-stag augmentation R_st[144] * sstag_ue[1]
 * Compare each piece against mfoil's totals to localize the 0.11 abs error.
 */

#include <H5Cpp.h>

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "scirocco/core/adjoint.h"
#include "scirocco/core/boundary_layer.h"
#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/forces.h"
#include "scirocco/core/wake.h"
#include "scirocco/properties/bl_properties.h"
#include "scirocco/scirocco.h"

using namespace scirocco;
using namespace Eigen;

static MatrixXd readMatrix(H5::H5File& f, const std::string& p) {
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
static VectorXd readVector(H5::H5File& f, const std::string& p) {
    H5::DataSet d = f.openDataSet(p);
    hsize_t dims[1];
    d.getSpace().getSimpleExtentDims(dims);
    std::vector<double> b(dims[0]);
    d.read(b.data(), H5::PredType::NATIVE_DOUBLE);
    return Map<VectorXd>(b.data(), dims[0]);
}

/**
 * @brief Entry point: decomposes R_V[144,195] into simi-station and xi-stag
 *        augmentation contributions and compares against mfoil reference.
 * @return 0 on success.
 */
int main() {
    const std::string path =
        std::string(SCIROCCO_TEST_DATA_DIR) + "/comprehensive_reference.h5";
    H5::H5File f(path, H5F_ACC_RDONLY);
    MatrixXd foil_x = readMatrix(f, "comprehensive/foil_x");
    MatrixXd U_ref = readMatrix(f, "comprehensive/U");
    VectorXd wgap_ref = readVector(f, "comprehensive/wgap");
    MatrixXd R_V_ref = readMatrix(f, "comprehensive/R_V");
    int N = foil_x.cols();

    Scirocco foil("0012", N - 1);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    auto& s = foil.GetState();
    s.foil.x = foil_x;
    s.foil.N = N;
    RecomputeArcLength(s);
    s.oper.viscous = true;
    s.oper.initbl = true;
    s.param.verb = 0;
    s.isol.AIC.resize(0, 0);
    s.isol.gamref.resize(0, 0);
    s.isol.gam.resize(0);
    SolveViscous(s);
    /* Override to the reference state and refresh sstag_ue */
    s.glob.U = U_ref;
    s.vsol.wgap = wgap_ref;
    {
        int I0 = s.isol.Istag[0], I1 = s.isol.Istag[1];
        double ue0 = s.glob.U(3, I0), ue1 = s.glob.U(3, I1);
        double S0 = s.foil.s(I0), S1 = s.foil.s(I1);
        double den = ue0 + ue1;
        s.isol.sstag_ue.resize(2);
        s.isol.sstag_ue(0) = ue1 * (S1 - S0) / (den * den);
        s.isol.sstag_ue(1) = -ue0 * (S1 - S0) / (den * den);
    }
    s.glob.realloc = true;
    BuildGlobSys(s);
    StoreAdjointFactorization(s);

    /* Now: decompose R_V[144, 195]. */
    int Istag0 = s.isol.Istag[0]; /* node 47 */
    int Istag1 = s.isol.Istag[1]; /* node 48 */
    std::cout << "Istag = [" << Istag0 << ", " << Istag1 << "]\n";
    std::cout << "sstag_ue = [" << s.isol.sstag_ue(0) << ", "
              << s.isol.sstag_ue(1) << "]\n";

    /* (1) Compute the simi-station R1_U at node Istag1=48 manually.
     * Replicate the BuildGlobSys simi-station path (coupled_solver.cpp 240+).
     */
    int si = 1; /* upper surface */
    const auto& Is = s.vsol.Is[si];
    Param param = BuildParam(s, si);
    VectorXd xi(Is.size());
    for (size_t k = 0; k < Is.size(); ++k) xi(k) = s.isol.xi(Is[k]);
    VectorXd Aux = VectorXd::Zero(Is.size());
    int ip0 = 0, ip1 = 1;
    Matrix42d Uin;
    Uin.col(0) = s.glob.U.col(Is[ip0]);
    Uin.col(1) = s.glob.U.col(Is[ip1]);
    Vector2d xin(xi(ip0), xi(ip1));
    auto stag = StagnationState(Uin, xin);
    param.turb = false;
    param.simi = true;
    Matrix42d Ust_pair;
    Ust_pair.col(0) = stag.Ust;
    Ust_pair.col(1) = stag.Ust;
    Vector2d xst_pair(stag.xst, stag.xst);
    Vector2d Aux_pair(Aux(ip0), Aux(ip0));
    auto res = ResidualStation(param, xst_pair, Ust_pair, Aux_pair);
    Matrix34d R1_Ust =
        res.lin_U.block<3, 4>(0, 0) + res.lin_U.block<3, 4>(0, 4);
    Matrix34d R1_U_node48 = R1_Ust * stag.Ust_U.block<4, 4>(0, 0);
    Matrix34d R1_U_node49 = R1_Ust * stag.Ust_U.block<4, 4>(0, 4);

    std::cout << "\n== Manual reconstruction of simi-station R1_U ==\n";
    std::cout << "R1_Ust (3x4):\n" << R1_Ust << "\n";
    std::cout << "stag.Ust:\n" << stag.Ust.transpose() << "\n";
    std::cout << "stag.Ust_U @ U[node48]:\n"
              << stag.Ust_U.block<4, 4>(0, 0) << "\n";
    std::cout << "R1_U_node48[0, 3] (R_mom wrt U[ue, node48]) = "
              << R1_U_node48(0, 3) << "\n";

    /* Print res.lin_U in detail for the simi station */
    std::cout.precision(15);
    std::cout << "\n== res.lin_U at simi station (3x8) ==\n";
    std::cout << res.lin_U << "\n";

    /* ∂Rmom/∂Ust[3] from both U1 and U2 contributions */
    std::cout << "R1_Ust(0, 3) = res.lin_U(0,3) + res.lin_U(0,7) = "
              << res.lin_U(0, 3) + res.lin_U(0, 7) << "\n";
    std::cout << "  res.lin_U(0,3) = " << res.lin_U(0, 3) << "\n";
    std::cout << "  res.lin_U(0,7) = " << res.lin_U(0, 7) << "\n";

    /* For comparison: probe what mfoil would have computed.
     * Trace: ∂Rmom/∂Ust[3] = (H_U[3] + Hw_U[3] - Ms_U[3]) * uelog
     *                       + (2+H+Hw-Ms)*uelog_U[3] - 0.5*xlog*cfxt_U[3]
     * At simi: thlog_U=0, uelog=1, uelog_U=0, xlog=1.
     * → ∂Rmom/∂Ust[3] = (H_U[3] + Hw_U[3] - Ms_U[3]) - 0.5 * cfxt_U[3]
     * H = ds/th, ∂H/∂ue = 0. H_U[3]=0.
     * Hw_U[3] = 0 (no ue dep).
     * Ms = M2, depends on uk and H0. For Ma=0: M2 = 0 → Ms_U[3] = 0.
     * So ∂Rmom/∂Ust[3] = -0.5 * cfxt_U[3].
     * cfxt = 0.25*cfxt1 + 0.5*cfxtm + 0.25*cfxt2.
     * For U1=U2=Ust, all three are GetCfxt(Ust, xst, param_simi_lam).
     * cfxt_U[3] (U1 part) = 0.25*cfxt1_U1[3] + 0.25*cfxtm_Um[3]*0.5*2  ...
     * Actually let me just print intermediate cf, cfxt values. */
    {
        Param p = param;
        p.turb = false;
        p.simi = true;
        auto [cf_at_Ust, cf_U_at_Ust] = GetCf(stag.Ust, p);
        auto [cfxt_at_Ust, cfxt_U_at_Ust, cfxt_x_at_Ust] =
            GetCfxt(stag.Ust, stag.xst, p);
        auto [Ret_at_Ust, Ret_U_at_Ust] = GetRet(stag.Ust, p);
        auto [Hk_at_Ust, Hk_U_at_Ust] = GetHk(stag.Ust, p);
        std::cout << "\n== At Ust=" << stag.Ust.transpose() << " ==\n";
        std::cout << "Ret  = " << Ret_at_Ust
                  << " (Ret_U = " << Ret_U_at_Ust.transpose() << ")\n";
        std::cout << "Hk   = " << Hk_at_Ust
                  << " (Hk_U  = " << Hk_U_at_Ust.transpose() << ")\n";
        std::cout << "cf   = " << cf_at_Ust
                  << " (cf_U  = " << cf_U_at_Ust.transpose() << ")\n";
        std::cout << "cfxt = " << cfxt_at_Ust
                  << " (cfxt_U=" << cfxt_U_at_Ust.transpose() << ")\n";
    }

    /* (2) xi-stag augmentation R_st * sstag_ue.
     * R_st[144] = (R_x * x_st)[144]. */
    int Nsys = s.glob.Nsys;
    int Nfoil = s.foil.N;
    VectorXd x_st(Nsys);
    for (int i = 0; i < Nfoil; ++i) x_st(i) = -s.isol.sgnue(i);
    for (int i = Nfoil; i < Nsys; ++i) x_st(i) = -1.0;
    VectorXd R_st = s.glob.R_x * x_st;
    std::cout << "\n== xi-stag augmentation ==\n";
    std::cout << "R_st[144] = " << R_st(144) << "\n";
    std::cout << "sstag_ue[1] (for Istag[1]=48) = " << s.isol.sstag_ue(1)
              << "\n";
    double aug = R_st(144) * s.isol.sstag_ue(1);
    std::cout << "augmentation = R_st[144] * sstag_ue[1] = " << aug << "\n";

    /* (3) Total scirocco R_V[144, 195]. */
    MatrixXd R_V_cpp = MatrixXd(s.glob.R_V);
    std::cout << "\n== Totals ==\n";
    std::cout << "scirocco  R_V[144, 195] = " << R_V_cpp(144, 195) << "\n";
    std::cout << "mfoil    R_V[144, 195] = " << R_V_ref(144, 195) << "\n";
    std::cout << "abs diff               = "
              << std::abs(R_V_cpp(144, 195) - R_V_ref(144, 195)) << "\n";
    std::cout << "manual R1_U[0,3] (no aug) = " << R1_U_node48(0, 3) << "\n";
    std::cout << "manual R1_U + aug         = " << R1_U_node48(0, 3) + aug
              << "\n";
    std::cout << "diff vs scirocco R_V       = "
              << std::abs(R1_U_node48(0, 3) + aug - R_V_cpp(144, 195)) << "\n";

    /* (4) Inferred mfoil R1_U(0, 3) = mfoil R_V[144, 195] - scirocco_aug.
     * (mfoil's aug should be same as scirocco's if sstag_ue and R_st agree) */
    double mfoil_R1U_inferred = R_V_ref(144, 195) - aug;
    std::cout
        << "\nInferred mfoil R1_U[0,3] = R_V_ref[144,195] - scirocco_aug = "
        << mfoil_R1U_inferred << "\n";
    std::cout << "Diff scirocco R1_U[0,3] vs inferred mfoil: "
              << std::abs(R1_U_node48(0, 3) - mfoil_R1U_inferred) << "  rel="
              << std::abs(R1_U_node48(0, 3) - mfoil_R1U_inferred) /
                     std::abs(mfoil_R1U_inferred)
              << "\n";

    return 0;
}
