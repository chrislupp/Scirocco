#!/usr/bin/env python3
"""mfoil-CS ∂J/∂Ma_partial at FIXED viscous U for ViscousMach config
(NACA 0012, α=3, Re=1e6, Ma=0.3). Holds U constant while perturbing Ma
(no Newton resolve), runs calc_force, takes Im(cl/cd/cm)/ε.

This is exactly the partial that cppfoil's analytic dJ_dMa block in
ViscousSensitivities Ma sub-block should produce.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'temp'))
import mfoil_cs  # noqa
from mfoil_cs import complex_step_mode  # noqa
import mfoil  # noqa

CS_STEP = 1e-30

M = mfoil.mfoil(naca='0012', npanel=100)
M.param.verb = 0
M.param.doplot = False
M.setoper(alpha=3.0, Re=1e6, Ma=0.3)
M.solve()
print(f'Base: cl={M.post.cl:.15g} cd={M.post.cd:.15g} cm={M.post.cm:.15g}')

# Now hold U fixed, perturb Ma, re-run init_thermo + calc_force.
# Promote all force-affecting arrays to complex.
M.glob.U = M.glob.U.astype(complex)
M.isol.gam = M.isol.gam.astype(complex) if hasattr(M.isol, 'gam') and M.isol.gam is not None else M.isol.gam
M.isol.uewi = M.isol.uewi.astype(complex) if hasattr(M.isol, 'uewi') and M.isol.uewi is not None else M.isol.uewi
M.isol.sgnue = M.isol.sgnue.astype(complex)
Ma0 = M.oper.Ma
with complex_step_mode():
    M.oper.Ma = Ma0 + 1j * CS_STEP
    mfoil.init_thermo(M)
    mfoil.calc_force(M)

dcl = float(np.imag(M.post.cl)) / CS_STEP
dcd = float(np.imag(M.post.cd)) / CS_STEP
dcm = float(np.imag(M.post.cm)) / CS_STEP
print(f'dJ/dMa at fixed U (CS):')
print(f'  dcl/dMa = {dcl:.15g}')
print(f'  dcd/dMa = {dcd:.15g}')
print(f'  dcm/dMa = {dcm:.15g}')
