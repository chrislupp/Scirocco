/**
 * @brief Compare scirocco's analytic adjoint gradient against scirocco's own
 * finite-difference gradient (no mfoil). If these agree at <1e-9 rel,
 * the analytic chain is internally consistent and the remaining floor
 * against the CS reference is a scirocco-mfoil primal divergence.
 * If they disagree, the analytic chain has a real bug.
 */

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "scirocco/core/adjoint.h"
#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/forces.h"
#include "scirocco/core/wake.h"
#include "scirocco/scirocco.h"

using namespace scirocco;
using namespace Eigen;

struct Force {
    double cl, cd, cm;
};

/**
 * @brief Solve the airfoil at the given conditions and return force coefficients.
 * @param alpha Angle of attack in degrees.
 * @param Ma Mach number.
 * @param Re Reynolds number.
 * @param foil_x_override Optional override for foil surface coordinates.
 * @return Force coefficients {cl, cd, cm}.
 */
static Force solveAt(double alpha, double Ma, double Re,
                     const MatrixXd* foil_x_override = nullptr) {
    Scirocco foil("0012", 100);
    foil.SetOper(alpha, -1.0, Re, Ma);
    auto& s = foil.GetState();
    if (foil_x_override) {
        s.foil.x = *foil_x_override;
        s.foil.N = foil_x_override->cols();
        RecomputeArcLength(s);
    }
    s.oper.viscous = true;
    s.oper.initbl = true;
    s.param.verb = 0;
    s.isol.AIC.resize(0, 0);
    s.isol.gamref.resize(0, 0);
    s.isol.gam.resize(0);
    SolveViscous(s);
    if (!s.glob.conv) {
        std::cerr << "Newton did not converge at alpha=" << alpha
                  << " Ma=" << Ma << " Re=" << Re << std::endl;
        return {NAN, NAN, NAN};
    }
    CalcForce(s);
    return {s.post.cl, s.post.cd, s.post.cm};
}

/**
 * @brief Entry point: compares analytic sensitivities against FD for alpha and
 *        shape perturbations.
 * @return 0 on success.
 */
int main() {
    /* -- Alpha FD vs analytic ---------------------------------------------- */
    std::cout << "=== Alpha FD vs analytic (Re=1e6, alpha=3) ===\n";
    {
        /* Analytic */
        Scirocco foil("0012", 100);
        foil.SetOper(3.0, -1.0, 1e6, 0.0);
        foil.Solve();
        auto sens = foil.ComputeSensitivities();
        double a_dcl = sens.dcl_dalpha;
        double a_dcd = sens.dcd_dalpha;
        double a_dcm = sens.dcm_dalpha;
        std::cout << "Analytic: dcl=" << a_dcl << " dcd=" << a_dcd
                  << " dcm=" << a_dcm << "\n";

        /* FD with several step sizes */
        for (double h : {1e-3, 1e-4, 1e-5, 1e-6}) {
            auto fp = solveAt(3.0 + h, 0.0, 1e6);
            auto fm = solveAt(3.0 - h, 0.0, 1e6);
            double fd_cl = (fp.cl - fm.cl) / (2.0 * h);
            double fd_cd = (fp.cd - fm.cd) / (2.0 * h);
            double fd_cm = (fp.cm - fm.cm) / (2.0 * h);
            std::cout << "h=" << h << ":  FD dcl=" << fd_cl
                      << " (rel=" << std::abs((a_dcl - fd_cl) / fd_cl) << ")"
                      << "  FD dcd=" << fd_cd
                      << " (rel=" << std::abs((a_dcd - fd_cd) / fd_cd) << ")"
                      << "  FD dcm=" << fd_cm
                      << " (rel=" << std::abs((a_dcm - fd_cm) / fd_cm) << ")\n";
        }
    }

    /* -- Shape FD at node 20, dim 0 (worst rel case) ------------------------ */
    std::cout << "\n=== Shape FD at node 20 dim 0 vs analytic (Re=1e6, "
                 "alpha=3) ===\n";
    {
        Scirocco foil("0012", 100);
        foil.SetOper(3.0, -1.0, 1e6, 0.0);
        auto& s = foil.GetState();
        s.oper.viscous = true;
        s.oper.initbl = true;
        s.param.verb = 0;
        s.isol.AIC.resize(0, 0);
        s.isol.gamref.resize(0, 0);
        s.isol.gam.resize(0);
        /* Snapshot foil.x BEFORE solve so FD perturbations have a clean
         * baseline */
        MatrixXd x_base = s.foil.x.leftCols(s.foil.N);
        /* For consistency with the sensitivity test, run RecomputeArcLength
         * (matches its convention) */
        RecomputeArcLength(s);
        SolveViscous(s);
        CalcForce(s);
        auto sens = ViscousSensitivities(s);
        double a_dcl = sens.dcl_dx(0, 20);
        double a_dcd = sens.dcd_dx(0, 20);
        double a_dcm = sens.dcm_dx(0, 20);
        std::cout << "Analytic node 20 dim 0: dcl=" << a_dcl << " dcd=" << a_dcd
                  << " dcm=" << a_dcm << "\n";

        for (double h : {1e-5, 1e-6, 1e-7, 1e-8}) {
            MatrixXd xp = x_base, xm = x_base;
            xp(0, 20) += h;
            xm(0, 20) -= h;
            auto fp = solveAt(3.0, 0.0, 1e6, &xp);
            auto fm = solveAt(3.0, 0.0, 1e6, &xm);
            double fd_cl = (fp.cl - fm.cl) / (2.0 * h);
            double fd_cd = (fp.cd - fm.cd) / (2.0 * h);
            double fd_cm = (fp.cm - fm.cm) / (2.0 * h);
            std::cout << "h=" << h << ":  FD dcl=" << fd_cl
                      << " (rel=" << std::abs((a_dcl - fd_cl) / fd_cl) << ")"
                      << "  FD dcd=" << fd_cd
                      << " (rel=" << std::abs((a_dcd - fd_cd) / fd_cd) << ")"
                      << "  FD dcm=" << fd_cm
                      << " (rel=" << std::abs((a_dcm - fd_cm) / fd_cm) << ")\n";
        }
    }
    return 0;
}
