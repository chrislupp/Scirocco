// Compute cppfoil's analytic dJ_dMa partial (at fixed U) and compare to
// mfoil-CS reference printed by extract_dJ_dMa_cs.py.

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/forces.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/cppfoil.h"
#include "cppfoil/properties/bl_properties.h"
#include "cppfoil/util/math_utils.h"

using namespace cppfoil;
using namespace Eigen;

// Mirror of ComputeCpDerivsMa from adjoint.cpp (anonymous static there).
static void ComputeCpDerivsMa(double ue, const Param& param, double Ma,
                              double& dcp_dMa, double& d2cp_duedMa) {
    double Vinf = std::max(param.Vinf, 1e-10);
    double cp_inc = 1.0 - (ue / Vinf) * (ue / Vinf);
    double dcp_due_inc = -2.0 * ue / (Vinf * Vinf);
    if (Ma <= 0 || param.KTl <= 0) {
        dcp_dMa = 0;
        d2cp_duedMa = 0;
        return;
    }
    double b = param.KTb;
    double l = param.KTl;
    double db = -Ma / b;
    double opb = 1.0 + b;
    double dl =
        2.0 * Ma / (opb * opb) + 2.0 * Ma * Ma * Ma / (opb * opb * opb * b);
    double den = b + 0.5 * l * (1.0 + b) * cp_inc;
    double den_cp = 0.5 * l * (1.0 + b);
    double cp = cp_inc / den;
    double dden_dMa = db + 0.5 * (dl * opb + l * db) * cp_inc;
    double dden_cp_dMa = 0.5 * (dl * opb + l * db);
    dcp_dMa = -cp / den * dden_dMa;
    double dfac_dMa = (-(dcp_dMa * den_cp + cp * dden_cp_dMa) * den -
                       (1.0 - cp * den_cp) * dden_dMa) /
                      (den * den);
    d2cp_duedMa = dcp_due_inc * dfac_dMa;
}

int main() {
    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, 0.3);
    auto& s = foil.GetState();
    s.oper.viscous = true;
    s.oper.initbl = true;
    s.param.verb = 0;
    SolveViscous(s);
    CalcForce(s);
    std::cout << "Base: cl=" << s.post.cl << " cd=" << s.post.cd
              << " cm=" << s.post.cm << "\n";

    // Replicate the dJ_dMa(0) = dcl/dMa code from ViscousSensitivities Ma block
    int N = s.foil.N;
    double Ma = s.oper.Ma;
    double chord = s.geom.chord;
    double alpha = s.oper.alpha;
    double cosa = Cosd(alpha);
    double sina = Sind(alpha);
    Vector2d xref = s.geom.xref;
    VectorXd ue = s.glob.U.row(3).head(N).transpose();
    double dcl = 0, dcm = 0, dcdpi = 0;
    for (int i0 = 1; i0 <= N; ++i0) {
        int i = (i0 == N) ? 0 : i0;
        int ip = (i0 == N) ? N - 1 : i0 - 1;
        Vector2d x1 = s.foil.x.col(ip);
        Vector2d x2 = s.foil.x.col(i);
        Vector2d dxv = x2 - x1;
        Vector2d dx1 = x1 - xref;
        Vector2d dx2 = x2 - xref;
        double dx_lift = -dxv(0) * cosa - dxv(1) * sina;
        double dz = dxv(1) * cosa - dxv(0) * sina;
        double dcp_ip, d2cp_ip, dcp_i, d2cp_i;
        ComputeCpDerivsMa(ue(ip), s.param, Ma, dcp_ip, d2cp_ip);
        ComputeCpDerivsMa(ue(i), s.param, Ma, dcp_i, d2cp_i);
        double dcpbar = 0.5 * (dcp_ip + dcp_i);
        dcl += dx_lift * dcpbar;
        dcdpi += dz * dcpbar;
        double dx1nds = dxv(0) * dx1(0) + dxv(1) * dx1(1);
        double dx2nds = dxv(0) * dx2(0) + dxv(1) * dx2(1);
        dcm += dcp_ip / 3.0 * dx1nds + dcp_ip / 6.0 * dx2nds +
               dcp_i / 6.0 * dx1nds + dcp_i / 3.0 * dx2nds;
    }
    std::cout << "Analytic dJ/dMa (cl/cm/cdpi panel + viscous Squire-Young for "
                 "cd):\n";
    std::cout << "  dcl/dMa     = " << dcl / chord << "\n";
    std::cout << "  dcm/dMa     = " << dcm / (chord * chord) << "\n";
    std::cout << "  dcdpi/dMa   = " << dcdpi / chord << "\n";

    // cd via Squire-Young
    if (!s.vsol.Is[2].empty()) {
        int iw = s.vsol.Is[2].back();
        Vector4d Uw = s.glob.U.col(iw);
        auto [H, _H_U] = GetH(Uw);
        (void)_H_U;
        auto [uk, _uk_U] = GetUk(Uw(3), s.param);
        (void)_uk_U;
        double Vinf = s.param.Vinf;
        double KTb = s.param.KTb;
        double KTl = s.param.KTl;
        double opb = 1.0 + KTb;
        double KTb_dMa_v = -Ma / KTb;
        double KTl_dMa_v = 2.0 * Ma / (opb * opb) -
                           2.0 * Ma * Ma / (opb * opb * opb) * KTb_dMa_v;
        double ue_v = Uw(3);
        double xKT = ue_v / Vinf;
        double den_uk = 1.0 - KTl * xKT * xKT;
        double uk_dMa =
            ue_v * (xKT * xKT - 1.0) / (den_uk * den_uk) * KTl_dMa_v;
        double q = uk / Vinf;
        double p = 0.5 * (5.0 + H);
        double qp = std::pow(q, p);
        double cd_dMa = 2.0 * Uw(0) * p * qp / q * (uk_dMa / Vinf);
        std::cout << "  dcd/dMa (SY) = " << cd_dMa << "\n";
    }
    return 0;
}
