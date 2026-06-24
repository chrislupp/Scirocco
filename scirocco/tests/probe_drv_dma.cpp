/**
 * @brief Diagnostic: analytic ∂R/∂Ma vs 4th-order FD on Ma at fixed U.
 */
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "scirocco/core/adjoint.h"
#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/scirocco.h"
using namespace scirocco;
using namespace Eigen;
/**
 * @brief Entry point: compares analytic dR/dMa against 4th-order FD and
 *        reports worst-case entries.
 * @return 0 on success.
 */
int main() {
    Scirocco foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, 0.3);
    auto& s = foil.GetState();
    s.oper.viscous = true;
    s.oper.initbl = true;
    s.param.verb = 0;
    SolveViscous(s);
    int Nsys = s.glob.Nsys;
    /* Force analytic (not the FD-fallback DRVDMa) */
    BuildGlobSys(s);
    VectorXd dR_an = VectorXd::Zero(4 * Nsys);
    dR_an.head(3 * Nsys) = s.glob.dR_dMa_times_Ma / s.oper.Ma;
    double Ma0 = s.oper.Ma;
    double h = std::abs(Ma0) * 1e-2;
    auto eval_RBL = [&](double Ma) {
        FoilState sp = s;
        sp.oper.Ma = Ma;
        InitThermo(sp);
        sp.glob.realloc = true;
        BuildGlobSys(sp);
        VectorXd R(4 * Nsys);
        R.head(3 * Nsys) = sp.glob.R;
        R.tail(Nsys).setZero();
        return R;
    };
    VectorXd Rp2 = eval_RBL(Ma0 + 2 * h);
    VectorXd Rp1 = eval_RBL(Ma0 + h);
    VectorXd Rm1 = eval_RBL(Ma0 - h);
    VectorXd Rm2 = eval_RBL(Ma0 - 2 * h);
    VectorXd dR_fd = (-Rp2 + 8 * Rp1 - 8 * Rm1 + Rm2) / (12 * h);
    double max_abs = 0;
    int wi = 0;
    for (int i = 0; i < 3 * Nsys; ++i) {
        double e = std::abs(dR_an(i) - dR_fd(i));
        if (e > max_abs) {
            max_abs = e;
            wi = i;
        }
    }
    std::cout << "Worst R_BL entry: row " << wi << " station " << wi / 3
              << " eq " << wi % 3 << " analytic=" << dR_an(wi)
              << " FD=" << dR_fd(wi) << " abs=" << max_abs << "\n";
    int n_diff = 0;
    for (int i = 0; i < 3 * Nsys; ++i) {
        if (std::abs(dR_an(i) - dR_fd(i)) > 1e-7 && std::abs(dR_fd(i)) > 1e-8) {
            ++n_diff;
            if (n_diff <= 5) {
                std::cout << "  station " << i / 3 << " eq " << i % 3
                          << ": an=" << dR_an(i) << " fd=" << dR_fd(i) << "\n";
            }
        }
    }
    std::cout << "Total entries with significant diff: " << n_diff << "\n";
    return 0;
}
