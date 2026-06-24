/**
 * @brief Print scirocco's sstag_ue, sgnue, Istag, and the augmentation R_st values
 * at the comprehensive_reference linearization point.
 */

#include <H5Cpp.h>

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "scirocco/core/adjoint.h"
#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/forces.h"
#include "scirocco/core/wake.h"
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
 * @brief Entry point: prints sstag_ue and stagnation quantities before and after
 *        overriding U with the mfoil reference state.
 * @return 0 on success.
 */
int main() {
    const std::string path =
        std::string(SCIROCCO_TEST_DATA_DIR) + "/comprehensive_reference.h5";
    H5::H5File f(path, H5F_ACC_RDONLY);
    MatrixXd foil_x = readMatrix(f, "comprehensive/foil_x");
    MatrixXd U_ref = readMatrix(f, "comprehensive/U");
    VectorXd wgap_ref = readVector(f, "comprehensive/wgap");
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

    std::cout << "After SolveViscous (scirocco's converged at cumulative):\n";
    std::cout << "  Istag = [" << s.isol.Istag[0] << ", " << s.isol.Istag[1]
              << "]\n";
    std::cout << "  sstag = " << s.isol.sstag << "\n";
    std::cout << "  sstag_ue = [" << s.isol.sstag_ue(0) << ", "
              << s.isol.sstag_ue(1) << "]\n";
    std::cout << "  sstag_g  = [" << s.isol.sstag_g(0) << ", "
              << s.isol.sstag_g(1) << "]\n";
    std::cout << "  ue at Istag = [" << s.glob.U(3, s.isol.Istag[0]) << ", "
              << s.glob.U(3, s.isol.Istag[1]) << "]\n";
    std::cout << "  gam at Istag = [" << s.isol.gam(s.isol.Istag[0]) << ", "
              << s.isol.gam(s.isol.Istag[1]) << "]\n";

    /* Compute sstag_ue manually from the formula using ue values */
    int I0 = s.isol.Istag[0];
    int I1 = s.isol.Istag[1];
    double ue0 = s.glob.U(3, I0);
    double ue1 = s.glob.U(3, I1);
    double S0 = s.foil.s(I0);
    double S1 = s.foil.s(I1);
    double den = ue0 + ue1;
    std::cout << "\nManual sstag_ue from converged ue, sstag formula:\n";
    std::cout << "  (ue=ue1*(S1-S0)/(den^2)):    "
              << ue1 * (S1 - S0) / (den * den) << "\n";
    std::cout << "  (ue=-ue0*(S1-S0)/(den^2)):   "
              << -ue0 * (S1 - S0) / (den * den) << "\n";

    /* sstag_g for comparison: it's the gamma-based equivalent */
    double g0 = s.isol.gam(I0);
    double g1 = s.isol.gam(I1);
    double gden = g1 - g0;
    std::cout << "\nManual sstag_g from gam, formula:\n";
    std::cout << "  g1*(S0-S1)/(gden^2):   " << g1 * (S0 - S1) / (gden * gden)
              << "\n";
    std::cout << "  -g0*(S0-S1)/(gden^2):  " << -g0 * (S0 - S1) / (gden * gden)
              << "\n";

    /* Now override U with mfoil's U_ref and recompute (mimicking the
     * residual_deriv_test pattern). Does sstag_ue stay the same? */
    s.glob.U = U_ref;
    s.vsol.wgap = wgap_ref;
    s.glob.realloc = true;
    BuildGlobSys(s);
    std::cout
        << "\nAfter override U=U_ref + BuildGlobSys (no StagPointMove):\n";
    std::cout << "  Istag = [" << s.isol.Istag[0] << ", " << s.isol.Istag[1]
              << "]\n";
    std::cout << "  sstag_ue = [" << s.isol.sstag_ue(0) << ", "
              << s.isol.sstag_ue(1) << "]\n";

    /* Manual recompute from U_ref */
    ue0 = s.glob.U(3, I0);
    ue1 = s.glob.U(3, I1);
    den = ue0 + ue1;
    std::cout << "\nManual sstag_ue from U_ref ue, sstag formula:\n";
    std::cout << "  (ue=ue1*(S1-S0)/(den^2)):    "
              << ue1 * (S1 - S0) / (den * den) << "\n";
    std::cout << "  (ue=-ue0*(S1-S0)/(den^2)):   "
              << -ue0 * (S1 - S0) / (den * den) << "\n";

    return 0;
}
