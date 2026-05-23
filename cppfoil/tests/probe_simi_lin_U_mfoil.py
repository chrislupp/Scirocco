#!/usr/bin/env python3
"""Probe mfoil's residual_station at the simi station with the same Ust
that cppfoil computes. Compare R1_Ust(0, 3) to cppfoil's value (228512)
vs the inferred mfoil value (228365).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'mfoil'))
import mfoil  # noqa: E402


def recompute_geom(M):
    x = M.foil.x
    N = M.foil.N
    ds = np.sqrt((x[0, 1:] - x[0, :-1])**2 + (x[1, 1:] - x[1, :-1])**2)
    s = np.zeros(N, dtype=x.dtype)
    s[1:] = np.cumsum(ds)
    M.foil.s = s


# Replicate the comprehensive_reference setup.
M = mfoil.mfoil(naca='0012', npanel=100)
M.param.verb = 0
M.param.doplot = False
M.setoper(alpha=3.0, Re=1e6, Ma=0.0)
M.solve()
recompute_geom(M)
mfoil.build_gamma(M, 3.0)
mfoil.build_wake(M)
mfoil.stagpoint_find(M)
mfoil.identify_surfaces(M)
mfoil.set_wake_gap(M)
mfoil.calc_ue_m(M)
mfoil.build_glob_sys(M)

# Upper surface = si=1
si = 1
Is = M.vsol.Is[si]
xi = M.isol.xi[Is]
# Aux for foil (si<2): zeros
Aux = np.zeros(len(Is))

param = mfoil.build_param(M, si)
i0 = 0

# Stagnation state from U[:, [0,1]] of upper surface
Uin = M.glob.U[:, [Is[i0], Is[i0+1]]]
xin = xi[[i0, i0+1]]
Ust, Ust_U, Ust_x, xst = mfoil.stagnation_state(Uin, xin)

print(f"Ust = {Ust}")
print(f"xst = {xst}")
print(f"Ust_U[:,:4] (deriv wrt U[node48]):\n{Ust_U[:,:4]}")

param.turb, param.simi = False, True
Ust_pair = np.stack((Ust, Ust), axis=-1)  # shape (4, 2)
xst_pair = np.r_[xst, xst]
Aux_pair = Aux[[i0, i0]]

R1, R1_Ut, R1_x = mfoil.residual_station(param, xst_pair, Ust_pair, Aux_pair)

print(f"\nres.value = {R1}")
print(f"res.lin_U (3x8):\n{R1_Ut}")
print(f"res.lin_U(0, 3) [Rmom wrt Ust1[3]] = {R1_Ut[0, 3]:.15g}")
print(f"res.lin_U(0, 7) [Rmom wrt Ust2[3]] = {R1_Ut[0, 7]:.15g}")
print(f"R1_Ust(0, 3) = sum = {R1_Ut[0, 3] + R1_Ut[0, 7]:.15g}")

R1_Ust = R1_Ut[:, :4] + R1_Ut[:, 4:]
R1_U_node48 = R1_Ust @ Ust_U[:, :4]
print(f"\nR1_U_node48[0, 3] = {R1_U_node48[0, 3]:.15g}")

# Also print intermediate closure values at Ust
cf, cf_U = mfoil.get_cf(Ust, param)
Ret, Ret_U = mfoil.get_Ret(Ust, param)
Hk, Hk_U = mfoil.get_Hk(Ust, param)
cfxt, cfxt_U, cfxt_x = mfoil.get_cfxt(Ust, xst, param)
print(f"\n== Closures at Ust ==")
print(f"Ret  = {Ret}")
print(f"Hk   = {Hk}")
print(f"cf   = {cf} (cf_U = {cf_U})")
print(f"cfxt = {cfxt} (cfxt_U = {cfxt_U})")

# Now dump augmentation values from mfoil
print(f"\n=== mfoil augmentation ===")
print(f"mfoil sstag_ue = {M.isol.sstag_ue}")
print(f"mfoil Istag    = {M.isol.Istag}")

Nsys = M.glob.Nsys
N = M.foil.N
x_st = -M.isol.sgnue
x_st = np.concatenate((x_st, -np.ones(M.wake.N)))
# R_x is sparse — use mfoil's current R_x (built by build_glob_sys above)
R_st = M.glob.R_x @ x_st
print(f"R_st[144] = {R_st[144]:.15g}")
print(f"aug = R_st[144] * sstag_ue[1] = {R_st[144] * M.isol.sstag_ue[1]:.15g}")
print(f"total R_V[144, 195] = R1_U + aug = "
      f"{171.591126581409 + R_st[144] * M.isol.sstag_ue[1]:.15g}")

