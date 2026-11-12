"""Complex-step alpha sweep: cut-cell (mfoil) vs captured transition.

Usage: python alpha_sweep.py MODE PARAMS A0 A1 N [NPANEL] > out.txt
  MODE   'cutcell' or 'captured'
  PARAMS JSON dict of M.param attributes, e.g. '{"capt_L": 80, "capt_b": 0}',
         plus optional "naca" (default "2412") and "Re" (default 1e6)
Prints alpha, cd, dcd/dalpha (complex step), cl, converged.
"""
import os, sys, warnings
os.environ['OMP_NUM_THREADS'] = '1'; os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'mfoil'))
import numpy as np
from multiprocessing import Pool
H = 1e-20
def point(args):
    a, mode, cfg, npanel = args
    import mfoil as mf, mfoil_cs, mfoil_captured as mc
    from mfoil_cs import complex_step_mode
    cfg = dict(cfg); naca = cfg.pop('naca', '2412'); Re = cfg.pop('Re', 1e6)
    M = mf.mfoil(naca=naca, npanel=npanel); M.param.verb = 0; M.param.doplot = False
    for k, v in cfg.items(): setattr(M.param, k, v)
    with complex_step_mode():
        M.setoper(alpha=a + 1j*H, Re=Re)
        M.solve()
        if mode == 'captured':
            mc.convert_state(M)
            for it in range(40):
                mc.build_glob_sys(M)
                if np.linalg.norm(np.real(M.glob.R)) < 1e-11 and np.linalg.norm(np.imag(M.glob.R)) < 1e-9*H: break
                mc.solve_glob(M); mc.update_state(M); mf.stagpoint_move(M)
            U5 = M.glob.U; M.glob.U = U5[0:4]; mf.calc_force(M)
    conv = bool(np.linalg.norm(np.real(M.glob.R)) < 1e-10) if mode == 'captured' else True
    return a, float(np.real(M.post.cd)), float(np.imag(M.post.cd))/H, float(np.real(M.post.cl)), conv
if __name__ == '__main__':
    import json
    mode = sys.argv[1]; cfg = json.loads(sys.argv[2]); a0, a1, n = float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
    npanel = int(sys.argv[6]) if len(sys.argv) > 6 else 160
    al = np.linspace(a0, a1, n+1)
    with Pool(10) as p:
        res = p.map(point, [(a, mode, cfg, npanel) for a in al])
    for a, cd, dcd, cl, conv in res:
        print('%.4f %.9f %.6e %.6f %d' % (a, cd, dcd, cl, conv))
