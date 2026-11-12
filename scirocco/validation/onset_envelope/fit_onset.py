"""Fit the envelope onset ramp (OnsetModel::Envelope) to the stability data.

Reads data/env_*.json from envelope.py. For each Falkner-Skan profile it fits
the envelope asymptote N = s (Re_theta - R0) over 3 < N < 9 and normalizes the
local envelope slope, r = (dN/dRe_theta) / s, against z = log10(Re_theta/R0).
The model is the gamma CDF rf = P(3, k (z - zc)), with k fixed by
E[10^(zc + T)] = 1 (T ~ Gamma(3, k)) so the asymptote is unchanged. A
per-profile zc is fitted, then zc = c0 + c1 / (H - 1).

Usage: python fit_onset.py [--docs path/to/onset.json]
"""
import glob, json, os, sys, warnings
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammainc

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
N_ORDER = 3


def k_of(zc):
    return np.log(10.0)/(1.0 - 10.0**(zc/N_ORDER))


def ramp(z, zc):
    return gammainc(N_ORDER, np.maximum(k_of(zc)*(z - zc), 0.0))


def cubic(z, dl=0.1):
    rn = np.clip((z + dl)/(2*dl), 0, 1)
    return 3*rn**2 - 2*rn**3


def correlation(H):
    """mfoil/XFOIL lrc(Hk) and asymptotic slope da = dN/dRe_theta."""
    Hmi = 1.0/(H - 1.0)
    lrc = 2.492*Hmi**0.43 + 0.7*(np.tanh(14*Hmi - 9.24) + 1)
    da = 0.028*(H - 1) - 0.0345*np.exp(-(3.87*Hmi - 2.52)**2)
    return lrc, da


def load():
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, 'data', 'env_*.json'))):
        d = json.load(open(f))
        R, N = np.array(d['Ret']), np.array(d['N_env'])
        sel = (N > 3) & (N < 9)
        if sel.sum() < 8:
            continue
        s, b = np.polyfit(R[sel], N[sel], 1)
        R0 = -b/s
        z = np.log10(R/R0)
        r = np.gradient(N, R)/s
        m = (z > -0.45) & (z < 0.45) & (N < 9.5)
        out.append(dict(beta=d['beta'], H=d['H'], R=R, N=N, s=s, R0=R0,
                        zcrit=np.log10(d['Re_theta_crit']/R0), z=z[m], r=r[m]))
    return sorted(out, key=lambda c: c['H'])


def main():
    cases = load()
    print(' beta     H    slope/da  log10R0-lrc  zcrit   zc_fit  rms_env  rms_cubic')
    for c in cases:
        o = minimize_scalar(lambda zc: np.mean((ramp(c['z'], zc) - c['r'])**2),
                            bounds=(-0.7, -0.05), method='bounded')
        c['zc'] = o.x
        lrc, da = correlation(c['H'])
        print('%+.3f  %.3f   %.3f      %+.3f      %+.3f  %+.3f   %.3f    %.3f' % (
            c['beta'], c['H'], c['s']/da, np.log10(c['R0']) - lrc, c['zcrit'],
            o.x, np.sqrt(o.fun), np.sqrt(np.mean((cubic(c['z']) - c['r'])**2))))
    Hmi = np.array([1/(c['H'] - 1) for c in cases])
    c1, c0 = np.polyfit(Hmi, [c['zc'] for c in cases], 1)
    res = [c['zc'] - (c0 + c1*h) for c, h in zip(cases, Hmi)]
    print('\nzc = %.4f + %.4f / (Hk - 1)   (max residual %.3f)' % (c0, c1, max(map(abs, res))))
    if '--docs' in sys.argv:
        path = sys.argv[sys.argv.index('--docs') + 1]
        zs = np.linspace(-0.4, 0.45, 86)
        H0 = 2.591  # Blasius: model curves at this shape factor
        zc0 = c0 + c1/(H0 - 1)
        doc = dict(z=zs.round(4).tolist(), cubic=cubic(zs).round(4).tolist(),
                   envelope=ramp(zs, zc0).round(4).tolist(), H_model=H0,
                   stability=[dict(H=round(c['H'], 3), z=c['z'][::3].round(4).tolist(),
                                   r=c['r'][::3].round(4).tolist())
                              for c in cases if c['beta'] in (0.2, 0.0, -0.14)])
        json.dump(doc, open(path, 'w'), separators=(',', ':'))
        print('wrote', path)


if __name__ == '__main__':
    main()
