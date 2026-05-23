// Compare cppfoil's computed adjoint vectors λ_cl, λ_cd to mfoil's
// CS-derived ones (stored in comprehensive_reference.h5). If λ matches
// at <1e-9 rel, the bug is in the VJP assembly. If λ differs at 1e-5,
// the bug is in R_V (either the assembly or the LU solve precision).

#include <H5Cpp.h>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "cppfoil/core/adjoint.h"
#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/forces.h"
#include "cppfoil/core/wake.h"
#include "cppfoil/cppfoil.h"

using namespace cppfoil;
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

int main() {
    const char* path = "/Users/clupp/tools/analysis/ProjectTrysil/cppfoil/"
                       "tests/comprehensive_reference.h5";
    H5::H5File f(path, H5F_ACC_RDONLY);
    MatrixXd foil_x = readMatrix(f, "comprehensive/foil_x");
    MatrixXd U_ref = readMatrix(f, "comprehensive/U");
    VectorXd wgap_ref = readVector(f, "comprehensive/wgap");
    VectorXd lambda_cl_ref = readVector(f, "comprehensive/lambda_cl");
    VectorXd lambda_cd_ref = readVector(f, "comprehensive/lambda_cd");

    int N = foil_x.cols();
    Cppfoil foil("0012", N - 1);
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
    CalcForce(s);
    CalcForceLin(s);

    // Override U/wgap to match HDF5. CRITICAL: also refresh sstag_ue from
    // the override U's ue values, otherwise the xi-stag chain in R_U uses
    // stale sstag_ue from cppfoil's own Newton convergence (which lives at
    // a different state than U_ref) and R_V differs from mfoil's at ~1e-5.
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

    // Re-populate the output partials at the over-ridden U (cl_V, cm_V, cd_V
    // are functions of U). Note: CalcForce was already called above; we
    // re-evaluate at U_ref to be consistent with the lambda solve.
    CalcForce(s);
    CalcForceLin(s);

    VectorXd lambda_cl = SolveViscousAdjoint(s, s.post.cl_V);
    VectorXd lambda_cd = SolveViscousAdjoint(s, s.post.cd_V);

    std::cout << "Nsys=" << s.glob.Nsys << " 4*Nsys=" << 4 * s.glob.Nsys
              << " lambda_cl_ref.size=" << lambda_cl_ref.size() << "\n";

    auto compare = [&](const VectorXd& a, const VectorXd& b,
                       const std::string& name) {
        if (a.size() != b.size()) {
            std::cout << name << ": SIZE MISMATCH " << a.size() << " vs "
                      << b.size() << "\n";
            return;
        }
        double max_abs = 0, max_rel = 0;
        int wi_abs = 0, wi_rel = 0;
        for (int i = 0; i < a.size(); ++i) {
            double e = std::abs(a(i) - b(i));
            if (e > max_abs) { max_abs = e; wi_abs = i; }
            double mag = std::max(std::abs(a(i)), std::abs(b(i)));
            if (mag > 1e-6) {
                double rel = e / mag;
                if (rel > max_rel) { max_rel = rel; wi_rel = i; }
            }
        }
        std::cout << name << ": max abs=" << max_abs << " at i=" << wi_abs
                  << " (cpp=" << a(wi_abs) << " ref=" << b(wi_abs) << ")\n";
        std::cout << name << ": max rel (|val|>1e-6) = " << max_rel
                  << " at i=" << wi_rel
                  << " (cpp=" << a(wi_rel) << " ref=" << b(wi_rel) << ")\n";
    };
    compare(lambda_cl, lambda_cl_ref, "lambda_cl");
    compare(lambda_cd, lambda_cd_ref, "lambda_cd");

    return 0;
}
