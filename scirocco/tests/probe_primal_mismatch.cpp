// Probe: compare scirocco's Newton-converged viscous state to mfoil's CS
// reference state. Diagnoses why end-to-end shape sensitivities are stuck
// at 1e-5 rel even though all chain partials match at 1e-12.

#include <H5Cpp.h>

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

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
static double readScalar(H5::H5File& f, const std::string& p) {
    H5::DataSet d = f.openDataSet(p);
    double v;
    d.read(&v, H5::PredType::NATIVE_DOUBLE);
    return v;
}

int main() {
    const std::string path = std::string(SCIROCCO_TEST_DATA_DIR) + "/comprehensive_reference.h5";
    H5::H5File f(path, H5F_ACC_RDONLY);
    MatrixXd foil_x = readMatrix(f, "comprehensive/foil_x");
    MatrixXd U_ref = readMatrix(f, "comprehensive/U");
    VectorXd wgap_ref = readVector(f, "comprehensive/wgap");
    VectorXd xi_ref = readVector(f, "comprehensive/xi");
    VectorXd gam_ref = readVector(f, "comprehensive/gam");
    double cl_ref = readScalar(f, "comprehensive/cl_base");
    double cd_ref = readScalar(f, "comprehensive/cd_base");
    double cm_ref = readScalar(f, "comprehensive/cm_base");

    int N = foil_x.cols();
    Scirocco foil("0012", N - 1);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    auto& s = foil.GetState();
    s.foil.x = foil_x;
    s.foil.N = N;
    // NB: skipping RecomputeArcLength — comparing with mfoil's spline foil.s
    s.oper.viscous = true;
    s.oper.initbl = true;
    s.param.verb = 0;
    s.isol.AIC.resize(0, 0);
    s.isol.gamref.resize(0, 0);
    s.isol.gam.resize(0);
    SolveViscous(s);
    CalcForce(s);

    std::cout << "\n==== PRIMAL COMPARISON ====\n";
    std::cout << "cl: cpp=" << s.post.cl << " mfoil=" << cl_ref
              << " rel=" << std::abs(s.post.cl - cl_ref) / std::abs(cl_ref)
              << "\n";
    std::cout << "cd: cpp=" << s.post.cd << " mfoil=" << cd_ref
              << " rel=" << std::abs(s.post.cd - cd_ref) / std::abs(cd_ref)
              << "\n";
    std::cout << "cm: cpp=" << s.post.cm << " mfoil=" << cm_ref
              << " rel=" << std::abs(s.post.cm - cm_ref) / std::abs(cm_ref)
              << "\n";

    std::cout << "\n==== INVISCID GAMMA ====\n";
    double gam_max = 0, gam_rel_max = 0;
    int gam_iworst = 0;
    for (int i = 0; i < N; ++i) {
        double diff = std::abs(s.isol.gam(i) - gam_ref(i));
        double rel = diff / std::max(std::abs(gam_ref(i)), 1e-12);
        if (rel > gam_rel_max) {
            gam_rel_max = rel;
            gam_iworst = i;
        }
        gam_max = std::max(gam_max, diff);
    }
    std::cout << "gam max abs=" << gam_max << " worst rel=" << gam_rel_max
              << " at i=" << gam_iworst << " (cpp=" << s.isol.gam(gam_iworst)
              << " mfoil=" << gam_ref(gam_iworst) << ")\n";

    std::cout << "\n==== xi (arc-length) ====\n";
    double xi_max = 0;
    int xi_iworst = 0;
    for (int i = 0; i < xi_ref.size(); ++i) {
        double diff = std::abs(s.isol.xi(i) - xi_ref(i));
        if (diff > xi_max) {
            xi_max = diff;
            xi_iworst = i;
        }
    }
    std::cout << "xi max abs=" << xi_max << " at i=" << xi_iworst
              << " (cpp=" << s.isol.xi(xi_iworst)
              << " mfoil=" << xi_ref(xi_iworst) << ")\n";

    std::cout << "\n==== wgap ====\n";
    double wgap_max = 0;
    int wgap_iworst = 0;
    for (int i = 0; i < wgap_ref.size(); ++i) {
        double diff = std::abs(s.vsol.wgap(i) - wgap_ref(i));
        if (diff > wgap_max) {
            wgap_max = diff;
            wgap_iworst = i;
        }
    }
    std::cout << "wgap max abs=" << wgap_max << " at i=" << wgap_iworst << "\n";

    std::cout << "\n==== Istag / Itran ====\n";
    std::cout << "cpp Istag = [" << s.isol.Istag[0] << "," << s.isol.Istag[1]
              << "]\n";
    int n_turb = 0;
    for (bool t : s.vsol.turb)
        if (t) ++n_turb;
    std::cout << "cpp num turbulent stations = " << n_turb << "/"
              << s.vsol.turb.size() << "\n";
    // Print first turbulent station index per surface
    std::cout << "Xt (transition xi, [lower, upper]): " << s.vsol.Xt(0, 0)
              << ", " << s.vsol.Xt(1, 0) << "\n";
    std::cout << "Xt (transition x , [lower, upper]): " << s.vsol.Xt(0, 1)
              << ", " << s.vsol.Xt(1, 1) << "\n";

    std::cout << "\n==== U comparison (4 x Nsys) ====\n";
    int Nsys = U_ref.cols();
    std::cout << "Nsys (cpp/mfoil): " << s.glob.U.cols() << " / " << Nsys
              << "\n";
    if (s.glob.U.cols() != Nsys) {
        std::cout << "DIM MISMATCH; skipping\n";
        return 1;
    }
    const char* names[4] = {"theta", "deltastar", "sa", "ue"};
    for (int k = 0; k < 4; ++k) {
        double mx = 0, rmx = 0;
        int iw = 0;
        for (int i = 0; i < Nsys; ++i) {
            double diff = std::abs(s.glob.U(k, i) - U_ref(k, i));
            double rel = diff / std::max(std::abs(U_ref(k, i)), 1e-12);
            if (rel > rmx) {
                rmx = rel;
                iw = i;
                mx = diff;
            }
        }
        std::cout << "  U[" << names[k] << "]: worst rel=" << rmx
                  << " at i=" << iw << " (cpp=" << s.glob.U(k, iw)
                  << " mfoil=" << U_ref(k, iw) << " abs=" << mx << ")\n";
    }

    std::cout << "\n==== R_BL at U_ref (per-station breakdown) ====\n";
    MatrixXd U_cpp = s.glob.U;
    s.glob.U = U_ref;
    s.vsol.wgap = wgap_ref;
    s.glob.realloc = true;
    BuildGlobSys(s);
    int Nstat = s.glob.R.size() / 3;
    std::cout << "Top 20 stations by ||R_3||_inf at U=U_ref:\n";
    std::vector<std::pair<double, int>> ranked;
    ranked.reserve(Nstat);
    for (int j = 0; j < Nstat; ++j) {
        double m = 0;
        for (int k = 0; k < 3; ++k)
            m = std::max(m, std::abs(s.glob.R(3 * j + k)));
        ranked.emplace_back(m, j);
    }
    std::sort(ranked.begin(), ranked.end(),
              [](auto a, auto b) { return a.first > b.first; });
    for (int i = 0; i < std::min<int>(20, Nstat); ++i) {
        int j = ranked[i].second;
        std::cout << "  station " << j
                  << " |Rmass|=" << std::abs(s.glob.R(3 * j))
                  << " |Rmom|=" << std::abs(s.glob.R(3 * j + 1))
                  << " |Rlag|=" << std::abs(s.glob.R(3 * j + 2));
        // Find which surface
        int sur = -1, idx = -1;
        for (int si = 0; si < 3; ++si) {
            const auto& Is = s.vsol.Is[si];
            for (size_t kk = 0; kk < Is.size(); ++kk) {
                if (Is[kk] == j) {
                    sur = si;
                    idx = (int)kk;
                    break;
                }
            }
            if (sur >= 0) break;
        }
        const char* surname[3] = {"lower", "upper", "wake"};
        if (sur >= 0)
            std::cout << "  (surf=" << surname[sur] << " idx=" << idx
                      << " turb=" << (s.vsol.turb[j] ? "T" : "L") << ")";
        std::cout << "\n";
    }
    std::cout << "Total Nsys=" << Nstat << ", num turb at U=U_ref state="
              << std::count(s.vsol.turb.begin(), s.vsol.turb.end(), true)
              << "\n";

    return 0;
}
