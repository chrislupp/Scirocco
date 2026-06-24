/**
 * @brief Targeted ping at wake station 1 — break down Rshape and Rlag analytic
 * contributions vs FD per closure. Specifically tests Hs_dMa, cDi_dMa,
 * cDixt_dMa via FD on the public closures.
 */

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "scirocco/core/bl_chain_internals.h"
#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/properties/bl_properties.h"
#include "scirocco/scirocco.h"

using namespace scirocco;
using namespace Eigen;

/**
 * @brief Entry point: pings each closure's FD ∂/∂Ma at wake station 1 and
 *        compares against the analytic _dMa values.
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

    int N = s.foil.N;
    int Nw = s.wake.N;
    double Ma0 = s.oper.Ma;
    double h = std::abs(Ma0) * 1e-5;
    Param base_param = s.param;

    /* Wake station 1 = vsol.Is[2][1]. Use U at that node. */
    int wnode = s.vsol.Is[2][1];
    Vector4d U = s.glob.U.col(wnode);
    double x = s.isol.xi(wnode);

    /* Build param matching wake station 1 in ResidualStation (wake=true,
     * turb=true after BuildParam + StationParam) */
    Param p = base_param;
    p.wake = true;
    p.turb = true;
    p.simi = false;

    auto eval_FD = [&](auto fn) {
        auto eval_at = [&](double Ma_pert) {
            s.oper.Ma = Ma_pert;
            InitThermo(s);
            Param pp = s.param;
            pp.wake = true;
            pp.turb = true;
            pp.simi = false;
            return fn(U, pp, x);
        };
        double f_p2 = eval_at(Ma0 + 2 * h);
        double f_p1 = eval_at(Ma0 + h);
        double f_m1 = eval_at(Ma0 - h);
        double f_m2 = eval_at(Ma0 - 2 * h);
        s.oper.Ma = Ma0;
        InitThermo(s);
        return (-f_p2 + 8 * f_p1 - 8 * f_m1 + f_m2) / (12 * h);
    };

    std::cout << "Wake station 1 (foil node " << wnode << ") at Ma=" << Ma0
              << "\n";
    std::cout << "Closure FD ∂/∂Ma values:\n";
    std::cout << "  GetUk(U[3])     = "
              << eval_FD([](const Vector4d& U, const Param& p, double) {
                     return GetUk(U(3), p).value;
                 })
              << "\n";
    std::cout << "  GetMach2(U)     = "
              << eval_FD([](const Vector4d& U, const Param& p, double) {
                     return GetMach2(U, p).value;
                 })
              << "\n";
    std::cout << "  GetHk(U)        = "
              << eval_FD([](const Vector4d& U, const Param& p, double) {
                     return GetHk(U, p).value;
                 })
              << "\n";
    std::cout << "  GetHs(U,turb=T) = "
              << eval_FD([](const Vector4d& U, const Param& p, double) {
                     return GetHs(U, p).value;
                 })
              << "\n";
    Param p_lam = p;
    p_lam.turb = false;
    std::cout << "  GetHs(U,turb=F) = "
              << eval_FD([](const Vector4d& U, const Param& p, double) {
                     Param pl = p;
                     pl.turb = false;
                     return GetHs(U, pl).value;
                 })
              << "\n";
    std::cout << "  GetRet(U)       = "
              << eval_FD([](const Vector4d& U, const Param& p, double) {
                     return GetRet(U, p).value;
                 })
              << "\n";
    std::cout << "  GetCDi(U)       = "
              << eval_FD([](const Vector4d& U, const Param& p, double) {
                     return GetCDi(U, p).value;
                 })
              << "\n";
    std::cout << "  GetCDixt(U,x)   = "
              << eval_FD([](const Vector4d& U, const Param& p, double x) {
                     return GetCDixt(U, x, p).value;
                 })
              << "\n";
    std::cout << "  GetUs(U)        = "
              << eval_FD([](const Vector4d& U, const Param& p, double) {
                     return GetUs(U, p).value;
                 })
              << "\n";
    std::cout << "  GetDe(U)        = "
              << eval_FD([](const Vector4d& U, const Param& p, double) {
                     return GetDe(U, p).value;
                 })
              << "\n";
    std::cout << "  GetCteq(U)      = "
              << eval_FD([](const Vector4d& U, const Param& p, double) {
                     return GetCteq(U, p).value;
                 })
              << "\n";
    std::cout << "  GetUk @ Vinf*(ue/Vinf=" << U(3) / p.Vinf << ")\n";

    /* Print closure VALUES at baseline */
    auto [Hk_v, _hu] = GetHk(U, p);
    auto [Hs_v, _Hsu] = GetHs(U, p);
    auto [Ret_v, _Ru] = GetRet(U, p);
    auto [cDi_v, _cDiu] = GetCDi(U, p);
    auto [cDixt_v, _cDixtU, _cDixt_x] = GetCDixt(U, x, p);
    Param p_l = p;
    p_l.turb = false;
    auto [Hs_l, _Hslu] = GetHs(U, p_l);
    std::cout << "\nValues at base: Hk=" << Hk_v << " Hs(turb)=" << Hs_v
              << " Hs(lam)=" << Hs_l << " Ret=" << Ret_v << " cDi=" << cDi_v
              << " cDixt=" << cDixt_v << "\n";

    /* Analytic _dMa values for direct comparison */
    detail::ThermoMaDeriv d = detail::ComputeThermoMaDeriv(p, p.Vinf, 0.0);
    std::cout << "\nAnalytic _dMa (use turb=T param):\n";
    std::cout << "  GetUk_dMa     = " << detail::GetUk_dMa(U(3), p, d) << "\n";
    std::cout << "  GetMach2_dMa  = " << detail::GetMach2_dMa(U, p, d) << "\n";
    std::cout << "  GetHk_dMa     = " << detail::GetHk_dMa(U, p, d) << "\n";
    std::cout << "  GetHs_dMa(T)  = " << detail::GetHs_dMa(U, p, d) << "\n";
    Param pl = p;
    pl.turb = false;
    std::cout << "  GetHs_dMa(F)  = " << detail::GetHs_dMa(U, pl, d) << "\n";
    std::cout << "  GetRet_dMa    = " << detail::GetRet_dMa(U, p, d) << "\n";
    std::cout << "  GetCDi_dMa    = " << detail::GetCDi_dMa(U, p, d) << "\n";
    std::cout << "  GetCDixt_dMa  = " << detail::GetCDixt_dMa(U, x, p, d)
              << "\n";
    std::cout << "  GetUs_dMa     = " << detail::GetUs_dMa(U, p, d) << "\n";

    return 0;
}
