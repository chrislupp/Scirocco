// Compare analytic DRVDRe (Ret-chain) against 4th-order central FD on Re
// at fixed U. Diagnostic for the partial Re chain — finds which residual
// rows the analytic chain misses.

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "cppfoil/core/adjoint.h"
#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/wake.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/cppfoil.h"

using namespace cppfoil;
using namespace Eigen;

int main() {
    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    auto& s = foil.GetState();
    s.oper.viscous = true;
    s.oper.initbl = true;
    s.param.verb = 0;
    SolveViscous(s);

    int Nsys = s.glob.Nsys;
    VectorXd dR_analytic = DRVDRe(s);

    // 4th-order central FD on Re at fixed U (rebuild thermo each time).
    double Re0 = s.oper.Re;
    double h = std::abs(Re0) * 1e-2;
    auto eval_RBL = [&](double Re) {
        FoilState sp = s;
        sp.oper.Re = Re;
        InitThermo(sp);
        sp.glob.realloc = true;
        BuildGlobSys(sp);
        VectorXd R(4 * Nsys);
        R.head(3 * Nsys) = sp.glob.R;
        R.tail(Nsys).setZero();  // ue rows have no Re dep
        return R;
    };
    VectorXd Rp2 = eval_RBL(Re0 + 2 * h);
    VectorXd Rp1 = eval_RBL(Re0 + h);
    VectorXd Rm1 = eval_RBL(Re0 - h);
    VectorXd Rm2 = eval_RBL(Re0 - 2 * h);
    VectorXd dR_fd = (-Rp2 + 8 * Rp1 - 8 * Rm1 + Rm2) / (12 * h);

    std::cout << "Analytic vs FD ∂R/∂Re comparison at Re=" << Re0 << "\n";
    // Per-station summary: are SOME stations getting populated? Print first 10
    // non-zero analytic stations and first 10 non-zero FD-but-zero-analytic.
    std::cout << "\nStations where analytic > 1e-15:\n";
    int n_nonzero_analytic = 0;
    for (int j = 0; j < Nsys && n_nonzero_analytic < 15; ++j) {
        bool any_nz = false;
        for (int k = 0; k < 3; ++k) {
            if (std::abs(dR_analytic(3 * j + k)) > 1e-15) any_nz = true;
        }
        if (any_nz) {
            ++n_nonzero_analytic;
            std::cout << "  station " << j
                      << " R_mom:  an=" << dR_analytic(3 * j + 0)
                      << " fd=" << dR_fd(3 * j + 0)
                      << " | R_sh: an=" << dR_analytic(3 * j + 1)
                      << " fd=" << dR_fd(3 * j + 1)
                      << " | R_lag: an=" << dR_analytic(3 * j + 2)
                      << " fd=" << dR_fd(3 * j + 2) << "\n";
        }
    }
    std::cout << "\nStations where FD>1e-12 but analytic ≈0:\n";
    int n_gap = 0;
    for (int j = 0; j < Nsys && n_gap < 15; ++j) {
        for (int k = 0; k < 3; ++k) {
            int i = 3 * j + k;
            if (std::abs(dR_fd(i)) > 1e-12 &&
                std::abs(dR_analytic(i)) < 1e-15) {
                ++n_gap;
                const char* eqn[3] = {"mom", "shape", "lag"};
                std::cout << "  station " << j << " R_" << eqn[k]
                          << ": FD=" << dR_fd(i) << "\n";
                break;
            }
        }
    }
    std::cout << "\n";

    double max_abs = 0;
    int wi = 0;
    for (int i = 0; i < 3 * Nsys; ++i) {
        double e = std::abs(dR_analytic(i) - dR_fd(i));
        if (e > max_abs) { max_abs = e; wi = i; }
    }
    std::cout << "Worst entry: row " << wi
              << " (station " << wi / 3 << " eq " << wi % 3 << ")"
              << " analytic=" << dR_analytic(wi)
              << " FD=" << dR_fd(wi)
              << " abs=" << max_abs
              << " rel=" << max_abs / std::max(std::abs(dR_fd(wi)), 1e-300)
              << "\n";

    // Worst by equation type (mom, shape, lag)
    const char* eq[3] = {"mom", "shape", "lag"};
    for (int k = 0; k < 3; ++k) {
        double m = 0;
        int wjk = 0;
        for (int j = 0; j < Nsys; ++j) {
            int i = 3 * j + k;
            double e = std::abs(dR_analytic(i) - dR_fd(i));
            if (e > m) { m = e; wjk = j; }
        }
        std::cout << "Worst R_" << eq[k] << ": station " << wjk
                  << " analytic=" << dR_analytic(3 * wjk + k)
                  << " FD=" << dR_fd(3 * wjk + k)
                  << " abs=" << m << "\n";
    }

    // Surface info
    std::cout << "Istag=[" << s.isol.Istag[0] << ", " << s.isol.Istag[1] << "]\n";
    std::cout << "N_foil=" << s.foil.N << " N_wake=" << s.wake.N << "\n";
    // Where is transition?
    int nturb_lo = 0, nturb_up = 0;
    for (auto i : s.vsol.Is[0]) if (s.vsol.turb[i]) nturb_lo++;
    for (auto i : s.vsol.Is[1]) if (s.vsol.turb[i]) nturb_up++;
    std::cout << "Lower surface turb stations: " << nturb_lo << "/"
              << s.vsol.Is[0].size() << "\n";
    std::cout << "Upper surface turb stations: " << nturb_up << "/"
              << s.vsol.Is[1].size() << "\n";

    return 0;
}
