// Per-station ResidualLin.lin_Ma analytic vs FD on Ma. Identifies WHICH
// station has the analytic chain bug. Once a station is known to be off,
// the per-closure suspects can be hunted by varying which closures
// contribute at that station (laminar/turb/simi branch).

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

#include "scirocco/core/boundary_layer.h"
#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/scirocco.h"
#include "scirocco/properties/bl_properties.h"

using namespace scirocco;
using namespace Eigen;

int main() {
    Scirocco foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, 0.3);
    auto& s = foil.GetState();
    s.oper.viscous = true;
    s.oper.initbl = true;
    s.param.verb = 0;
    SolveViscous(s);
    int N = s.foil.N;
    int Nw = s.wake.N;
    int Nsys = N + Nw;
    double Ma0 = s.oper.Ma;
    double h = std::abs(Ma0) * 1e-4;

    // For each non-first station i = 1..N-1 on each surface, compute
    // ResidualStation analytic lin_Ma AND FD lin_Ma. The first station
    // (i=0) uses the similarity-station path; skip for now.

    std::cout << "Per-station lin_Ma vs FD on ResidualStation:\n";
    std::cout << "  (analytic via ResidualLin.lin_Ma; FD via "
                 "(R(Ma+h)-R(Ma-h))/(2h)·Ma at fixed U,x,Aux)\n\n";

    double worst_rel = 0;
    int worst_si = 0, worst_i = 0, worst_eq = 0;
    Vector3d worst_an, worst_fd;

    for (int si = 0; si < 3; ++si) {
        const auto& Is = s.vsol.Is[si];
        if (Is.empty()) continue;
        int Ns = static_cast<int>(Is.size());
        Param param = BuildParam(s, si);
        // Skip first point (simi or wake-init); just march
        for (int i = 1; i < Ns; ++i) {
            int ip = i - 1;
            // Skip transition stations (different residual call)
            bool tran = (Is[ip] < static_cast<int>(s.vsol.turb.size())) &&
                        (Is[i] < static_cast<int>(s.vsol.turb.size())) &&
                        (s.vsol.turb[Is[ip]] != s.vsol.turb[Is[i]]);
            if (tran) continue;
            StationParam(s, param, Is[i]);

            Matrix42d U;
            U.col(0) = s.glob.U.col(Is[ip]);
            U.col(1) = s.glob.U.col(Is[i]);
            Vector2d x(s.isol.xi(Is[ip]), s.isol.xi(Is[i]));
            Vector2d Aux = Vector2d::Zero();
            if (si == 2 && static_cast<int>(s.vsol.wgap.size()) > 0) {
                Aux(0) = s.vsol.wgap(ip);
                Aux(1) = s.vsol.wgap(i);
            }

            // Analytic from ResidualLin
            auto res = ResidualStation(param, x, U, Aux);
            Vector3d lin_Ma_an = res.lin_Ma;  // already Vector3d·Ma

            // FD: perturb Ma + InitThermo on a copy of state's Param,
            // then call ResidualStation with the perturbed param.
            // For the FD we need to update param.{KTb, KTl, H0, mu0, rho0,
            // Minf, gam, cps, muinf} since those are what closures read.
            auto eval_R = [&](double Ma) {
                FoilState sp = s;
                sp.oper.Ma = Ma;
                InitThermo(sp);
                Param p_p = BuildParam(sp, si);
                StationParam(sp, p_p, Is[i]);
                return ResidualStation(p_p, x, U, Aux).value;
            };
            Vector3d Rp2 = eval_R(Ma0 + 2 * h);
            Vector3d Rp1 = eval_R(Ma0 + h);
            Vector3d Rm1 = eval_R(Ma0 - h);
            Vector3d Rm2 = eval_R(Ma0 - 2 * h);
            Vector3d dR_dMa_fd = (-Rp2 + 8 * Rp1 - 8 * Rm1 + Rm2) / (12 * h);
            Vector3d lin_Ma_fd = dR_dMa_fd * Ma0;

            // Compare
            for (int eq = 0; eq < 3; ++eq) {
                double an = lin_Ma_an(eq);
                double fd = lin_Ma_fd(eq);
                double mag = std::max(std::abs(an), std::abs(fd));
                if (mag < 1e-10) continue;
                double rel = std::abs(an - fd) / mag;
                if (rel > worst_rel) {
                    worst_rel = rel;
                    worst_si = si;
                    worst_i = i;
                    worst_eq = eq;
                    worst_an = lin_Ma_an;
                    worst_fd = lin_Ma_fd;
                }
            }
        }
    }

    const char* sur[] = {"lower", "upper", "wake"};
    const char* eqn[] = {"mom  ", "shape", "lag  "};
    std::cout << "Worst rel error: surface=" << sur[worst_si]
              << " i=" << worst_i << " eq=" << eqn[worst_eq]
              << " rel=" << worst_rel << "\n";
    std::cout << "  analytic: [" << worst_an.transpose() << "]\n";
    std::cout << "  FD:       [" << worst_fd.transpose() << "]\n\n";

    // Per-station table: lin_Ma analytic vs FD for first 20 stations of
    // each surface
    std::cout << "Sample stations:\n";
    for (int si = 0; si < 3; ++si) {
        const auto& Is = s.vsol.Is[si];
        if (Is.empty()) continue;
        int Ns = static_cast<int>(Is.size());
        std::cout << "Surface " << sur[si] << " (" << Ns << " stations):\n";
        Param param = BuildParam(s, si);
        for (int i = 1; i < std::min(Ns, 10); ++i) {
            int ip = i - 1;
            bool tran = (Is[ip] < static_cast<int>(s.vsol.turb.size())) &&
                        (Is[i] < static_cast<int>(s.vsol.turb.size())) &&
                        (s.vsol.turb[Is[ip]] != s.vsol.turb[Is[i]]);
            if (tran) {
                std::cout << "  i=" << i << " [TRANSITION]\n";
                continue;
            }
            StationParam(s, param, Is[i]);
            Matrix42d U;
            U.col(0) = s.glob.U.col(Is[ip]);
            U.col(1) = s.glob.U.col(Is[i]);
            Vector2d x(s.isol.xi(Is[ip]), s.isol.xi(Is[i]));
            Vector2d Aux = Vector2d::Zero();
            if (si == 2 && static_cast<int>(s.vsol.wgap.size()) > 0) {
                Aux(0) = s.vsol.wgap(ip);
                Aux(1) = s.vsol.wgap(i);
            }
            auto res = ResidualStation(param, x, U, Aux);
            Vector3d an = res.lin_Ma;
            auto eval_R = [&](double Ma) {
                FoilState sp = s;
                sp.oper.Ma = Ma;
                InitThermo(sp);
                Param p_p = BuildParam(sp, si);
                StationParam(sp, p_p, Is[i]);
                return ResidualStation(p_p, x, U, Aux).value;
            };
            Vector3d Rp2 = eval_R(Ma0 + 2 * h);
            Vector3d Rp1 = eval_R(Ma0 + h);
            Vector3d Rm1 = eval_R(Ma0 - h);
            Vector3d Rm2 = eval_R(Ma0 - 2 * h);
            Vector3d fd = (-Rp2 + 8 * Rp1 - 8 * Rm1 + Rm2) / (12 * h) * Ma0;
            std::cout << "  i=" << i << " turb=" << (param.turb ? "T" : "L")
                      << " an=[" << an.transpose() << "]"
                      << " fd=[" << fd.transpose() << "]\n";
        }
    }
    return 0;
}
