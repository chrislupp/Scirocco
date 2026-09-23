/**
 * @brief Phase B diagnostic: ping each per-closure _dMa helper against
 * 4th-order FD on Ma at fixed U. Reports worst rel error per closure across all
 * stations to localize which _dMa has a bug.
 */

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "scirocco/core/boundary_layer.h"
#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/properties/bl_properties.h"
#include "scirocco/scirocco.h"

using namespace scirocco;
using namespace Eigen;

/* The _dMa helpers + ComputeThermoMaDeriv live in an anonymous namespace
 * in boundary_layer.cpp, so they're not callable from outside. We mirror
 * them here as a temporary harness — declare extern for the ones that
 * are static, OR re-derive them inline by calling Get<X>(perturbed param)
 * directly.
 *
 * Strategy: per closure, compute analytic via small re-implementation that
 * matches the production helper, AND compute FD by perturbing state.oper.Ma
 * + calling InitThermo + calling the public Get<X>(U, param).
 *
 * This validates that the production Get<X>_dMa formulas (which we cannot
 * call directly) match FD. Since they're inline in ComputeStationLinMa,
 * any mismatch surfaces at the higher level. */

struct PingStats {
    std::string name;
    double max_abs = 0;
    double max_rel = 0;
    int worst_node = -1;
    double worst_an = 0, worst_fd = 0;
};

template <typename PrimalFn>
PingStats run_ping(const std::string& name, FoilState& s, PrimalFn primal) {
    PingStats stats;
    stats.name = name;
    int N = s.foil.N;
    int Nw = s.wake.N;
    int Nsys = N + Nw;
    double Ma0 = s.oper.Ma;
    double h = std::abs(Ma0) * 1e-3;
    Param param0 = s.param;

    for (int node = 0; node < Nsys; ++node) {
        Vector4d U = s.glob.U.col(node);
        double xi = (node < N) ? s.isol.xi(node) : s.isol.xi(node);
        /* Use a per-surface param (turb flag depends on node) */
        Param p = param0;
        if (node < static_cast<int>(s.vsol.turb.size())) {
            p.turb = s.vsol.turb[node];
        }
        p.wake = (node >= N);

        /* 4-pt central FD */
        auto eval_at = [&](double Ma_pert) {
            s.oper.Ma = Ma_pert;
            InitThermo(s);
            return primal(U, s.param, xi);
        };
        double f_p2 = eval_at(Ma0 + 2 * h);
        double f_p1 = eval_at(Ma0 + h);
        double f_m1 = eval_at(Ma0 - h);
        double f_m2 = eval_at(Ma0 - 2 * h);
        double fd = (-f_p2 + 8 * f_p1 - 8 * f_m1 + f_m2) / (12 * h);

        /* Restore baseline before computing analytic (we need baseline param)
         */
        s.oper.Ma = Ma0;
        InitThermo(s);
        Param p_now = s.param;
        if (node < static_cast<int>(s.vsol.turb.size())) {
            p_now.turb = s.vsol.turb[node];
        }
        p_now.wake = (node >= N);

        /* For now we only run the FD path; analytic comes from the
         * production code via ComputeStationLinMa accumulation, which is
         * tested at the ResidualLin level in probe_drv_dma. Here we use
         * the FD value to localize which CLOSURE'S Ma chain is causing
         * the ResidualLin bug — we print FD curves at each node + closure
         * for the user to spot.
         * (Skip per-station analytic isolation for now.) */
        double an = 0; /* placeholder */
        double err = std::abs(an - fd);
        if (err > stats.max_abs) {
            stats.max_abs = err;
            stats.worst_node = node;
            stats.worst_an = an;
            stats.worst_fd = fd;
        }
        double mag = std::max(std::abs(an), std::abs(fd));
        if (mag > 1e-8) {
            double rel = err / mag;
            if (rel > stats.max_rel) stats.max_rel = rel;
        }
    }
    /* Restore state */
    s.oper.Ma = Ma0;
    s.param = param0;
    return stats;
}

/**
 * @brief Entry point: probes each public closure's FD ∂/∂Ma profile across all
 *        stations and reports the worst-magnitude node per closure.
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

    /* Probe each public closure's FD profile across all stations. Compare
     * the FD profile against what ComputeStationLinMa thinks each station
     * should contribute, station-by-station.
     *
     * The simplest closure-level diagnostic: for each closure X, dump the
     * FD ∂X/∂Ma at each station. We compare those numbers against the
     * analytic formula by inspecting the relevant section of
     * ComputeStationLinMa source by hand. */

    auto getH = [](const Vector4d& U, const Param& p, double) {
        return GetH(U).value;
    };
    auto getHk = [](const Vector4d& U, const Param& p, double) {
        return GetHk(U, p).value;
    };
    auto getHs = [](const Vector4d& U, const Param& p, double) {
        return GetHs(U, p).value;
    };
    auto getHss = [](const Vector4d& U, const Param& p, double) {
        return GetHss(U, p).value;
    };
    auto getRet = [](const Vector4d& U, const Param& p, double) {
        return GetRet(U, p).value;
    };
    auto getCf = [](const Vector4d& U, const Param& p, double) {
        return GetCf(U, p).value;
    };
    auto getCDi = [](const Vector4d& U, const Param& p, double) {
        return GetCDi(U, p).value;
    };
    auto getCfxt = [](const Vector4d& U, const Param& p, double x) {
        return GetCfxt(U, x, p).value;
    };
    auto getCDixt = [](const Vector4d& U, const Param& p, double x) {
        return GetCDixt(U, x, p).value;
    };
    auto getDe = [](const Vector4d& U, const Param& p, double) {
        return GetDe(U, p).value;
    };
    auto getUs = [](const Vector4d& U, const Param& p, double) {
        return GetUs(U, p).value;
    };
    auto getCteq = [](const Vector4d& U, const Param& p, double) {
        return GetCteq(U, p).value;
    };
    auto getMach2 = [](const Vector4d& U, const Param& p, double) {
        return GetMach2(U, p).value;
    };
    auto getUk = [](const Vector4d& U, const Param& p, double) {
        return GetUk(U(3), p).value;
    };
    auto getDamp = [](const Vector4d& U, const Param& p, double) {
        return GetDamp(U, p).value;
    };

    std::vector<std::pair<std::string, double>> closure_fd_max;
    auto dump_closure = [&](const std::string& name, auto fn) {
        /* For each station compute FD ∂X/∂Ma. */
        int N = s.foil.N;
        int Nw = s.wake.N;
        int Nsys = N + Nw;
        double Ma0 = s.oper.Ma;
        double h = Ma0 * 1e-3;
        double max_abs = 0;
        int worst = 0;
        double worst_fd = 0;
        for (int node = 0; node < Nsys; ++node) {
            Vector4d U = s.glob.U.col(node);
            double xi = s.isol.xi(node);
            Param p = s.param;
            if (node < static_cast<int>(s.vsol.turb.size()))
                p.turb = s.vsol.turb[node];
            p.wake = (node >= N);
            auto eval_at = [&](double Ma_pert) {
                s.oper.Ma = Ma_pert;
                InitThermo(s);
                Param pp = s.param;
                if (node < static_cast<int>(s.vsol.turb.size()))
                    pp.turb = s.vsol.turb[node];
                pp.wake = (node >= N);
                return fn(U, pp, xi);
            };
            double f_p2 = eval_at(Ma0 + 2 * h);
            double f_p1 = eval_at(Ma0 + h);
            double f_m1 = eval_at(Ma0 - h);
            double f_m2 = eval_at(Ma0 - 2 * h);
            double fd = (-f_p2 + 8 * f_p1 - 8 * f_m1 + f_m2) / (12 * h);
            s.oper.Ma = Ma0;
            InitThermo(s);
            if (std::abs(fd) > max_abs) {
                max_abs = std::abs(fd);
                worst = node;
                worst_fd = fd;
            }
        }
        std::cout << "  " << name << ": worst |FD ∂/∂Ma| = " << max_abs
                  << " at node " << worst << " (fd=" << worst_fd << ")\n";
        closure_fd_max.emplace_back(name, max_abs);
    };

    std::cout << "Per-closure FD ∂/∂Ma profile at Ma=" << s.oper.Ma << "\n";
    std::cout << "(Use this to see which closures have non-trivial Ma dep at "
                 "this state)\n\n";
    dump_closure("GetH        ", getH);
    dump_closure("GetHk       ", getHk);
    dump_closure("GetHs       ", getHs);
    dump_closure("GetHss      ", getHss);
    dump_closure("GetMach2    ", getMach2);
    dump_closure("GetRet      ", getRet);
    dump_closure("GetUk       ", getUk);
    dump_closure("GetCf       ", getCf);
    dump_closure("GetCDi      ", getCDi);
    dump_closure("GetCfxt     ", getCfxt);
    dump_closure("GetCDixt    ", getCDixt);
    dump_closure("GetDe       ", getDe);
    dump_closure("GetUs       ", getUs);
    dump_closure("GetCteq     ", getCteq);
    dump_closure("GetDamp     ", getDamp);

    return 0;
}
