"""Envelope N-factor vs Re_theta for one Falkner-Skan profile.

Locally parallel spatial Orr-Sommerfeld analysis. At each Re_dstar the
unstable band is found by anchoring on the temporally most unstable TS mode
and sweeping omega outward to both neutral points. Each dimensional
frequency is followed at constant F = omega_d nu / U_e^2, so omega = F Re_dstar,
and dN/dRe_theta = -alpha_i / (H theta_eta^2) for Falkner-Skan similarity.
"""
import os, sys, json, time
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS', '1')
import numpy as np
from multiprocessing import Pool
from scipy.optimize import brentq, minimize_scalar
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import orr_sommerfeld as m  # noqa: E402
S = None
def init(beta):
    global S; S = m.OS(beta, N=60, ymax=30, yi=2.5)

def band(R):
    """Return arrays (omega, -alpha_i) across the unstable band at Re_dstar=R."""
    r = minimize_scalar(lambda al: -S.ts_omega(al, R).imag, bounds=(0.03, 0.7), method='bounded', options={'xatol': 1e-6})
    w0 = S.ts_omega(r.x, R)
    if w0.imag <= 0:
        return [], []
    a0 = S.spatial_alpha(w0.real, R, r.x + 0j, w0)
    pts = {w0.real: -a0.imag}
    for sgn in (+1, -1):
        a = a0; w = w0.real; dw = 0.04*w0.real
        while True:
            wn = w + sgn*dw
            if wn <= 0: break
            try:
                an = S.spatial_alpha(wn, R, a, None if False else wn + 0j)
            except Exception:
                if dw < 1e-4*w0.real: break
                dw *= 0.5; continue
            if abs(an - a) > 0.2*abs(a):
                if dw < 1e-4*w0.real: break
                dw *= 0.5; continue
            a, w = an, wn
            pts[w] = -a.imag
            if -a.imag < 0: break
    ws = np.array(sorted(pts)); g = np.array([pts[x] for x in ws])
    return ws.tolist(), g.tolist()

def crit(beta):
    s = m.OS(beta, N=60, ymax=30, yi=2.5)
    def wi_max(Re):
        r = minimize_scalar(lambda al: -s.ts_omega(al, Re).imag, bounds=(0.03, 0.7), method='bounded', options={'xatol': 1e-7})
        return -r.fun
    Rs = np.geomspace(20.0, 60000.0, 70); v = [wi_max(R) for R in Rs]
    k = next(i for i in range(len(v)-1) if v[i] < 0 <= v[i+1])
    return brentq(wi_max, Rs[k], Rs[k+1], xtol=1e-6), s.H, s.theta_eta

if __name__ == '__main__':
    beta = float(sys.argv[1]); kmax = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
    t0 = time.time()
    Rcd, H, th = crit(beta)
    Ret = (Rcd/H)*np.geomspace(1.0005, kmax, 200)
    Red = H*Ret
    with Pool(10, initializer=init, initargs=(beta,)) as p:
        bands = p.map(band, Red.tolist())
    # F grid spanning all band frequencies
    Fmin = min(min(w)/R for (w, g), R in zip(bands, Red) if w)
    Fmax = max(max(w)/R for (w, g), R in zip(bands, Red) if w)
    Fs = np.geomspace(Fmin, Fmax, 400)
    G = np.zeros((len(Fs), len(Red)))
    for j, ((w, g), R) in enumerate(zip(bands, Red)):
        if not w: continue
        w = np.array(w); g = np.array(g)
        om = Fs*R
        inside = (om >= w[0]) & (om <= w[-1])
        G[inside, j] = np.maximum(np.interp(om[inside], w, g), 0.0)
    dNdRt = G/(H*th**2)
    N = np.zeros_like(dNdRt)
    N[:, 1:] = np.cumsum(0.5*(dNdRt[:, 1:] + dNdRt[:, :-1])*np.diff(Ret)[None, :], axis=1)
    env = N.max(axis=0)
    out = dict(beta=beta, H=H, theta_eta=th, Re_dstar_crit=Rcd, Re_theta_crit=Rcd/H,
               Ret=Ret.tolist(), N_env=env.tolist(),
               max_growth=[max(g) if g else 0.0 for (w, g) in bands])
    json.dump(out, open(os.path.join(HERE, 'data', 'env_%+.3f.json' % beta), 'w'))
    print('beta %+.3f H %.4f Re_d* crit %.1f Re_th crit %.1f  N_env(end) %.2f  %.0fs' % (
        beta, H, Rcd, Rcd/H, env[-1], time.time()-t0), flush=True)
