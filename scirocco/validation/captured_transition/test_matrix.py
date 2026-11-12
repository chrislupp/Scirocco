"""Cut cell vs captured transition over airfoils, Reynolds numbers, angles.

Usage: python test_matrix.py > data/matrix.txt
Columns: airfoil Re alpha | cut cl cd xt_up | captured-L60 conv cl cd |
captured-Narasimha conv cl cd zone10-90(upper, x/c)
"""
import os, sys, warnings, itertools
os.environ['OMP_NUM_THREADS'] = '1'; os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'mfoil'))
import numpy as np
from multiprocessing import Pool


def solve(naca, Re, alpha, mode):
    import mfoil as mf, mfoil_captured as mc
    M = mf.mfoil(naca=naca, npanel=160); M.param.verb = 0; M.param.doplot = False
    M.setoper(alpha=alpha, Re=Re)
    if mode == 'narasimha':
        M.param.capt_Lmodel = 'narasimha'
    try:
        if mode == 'cut':
            M.solve()
            up = M.vsol.Is[1]
            xt = next((M.foil.x[0, j] for j in up if M.vsol.turb[j]), np.nan)
            return (bool(M.glob.conv), M.post.cl, M.post.cd, xt, np.nan)
        conv = mc.solve_captured(M)
        g = mc.gamma(M.glob.U[mc.AMP], M.param)[0]
        up = [j for j in M.vsol.Is[1] if j < M.foil.N]
        z = [M.foil.x[0, j] for j in up if 0.1 < g[j] < 0.9]
        return (bool(conv), M.post.cl, M.post.cd,
                z[0] if z else np.nan, z[-1] if z else np.nan)
    except Exception:
        return (False, np.nan, np.nan, np.nan, np.nan)


def case(args):
    naca, Re, alpha = args
    return args, solve(naca, Re, alpha, 'cut'), solve(naca, Re, alpha, 'L60'), \
        solve(naca, Re, alpha, 'narasimha')


if __name__ == '__main__':
    cases = list(itertools.product(['0012', '2412', '4412'], [2e5, 1e6, 3e6], [0.0, 4.0, 8.0]))
    with Pool(10) as p:
        res = p.map(case, cases)
    print('# airfoil Re alpha | cut: conv cl cd xt | L60: conv cl cd dcd% | narasimha: conv cl cd dcd% zone')
    for (naca, Re, a), c, l, n in res:
        rel = lambda r: 100*(r[2] - c[2])/c[2]
        print('%s %.0e %4.1f | %d %.4f %.5f %.3f | %d %.4f %.5f %+5.1f | %d %.4f %.5f %+5.1f %.3f-%.3f' % (
            naca, Re, a, c[0], c[1], c[2], c[3], l[0], l[1], l[2], rel(l),
            n[0], n[1], n[2], rel(n), n[3], n[4]))
