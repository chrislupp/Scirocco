"""Calibrate the captured transition length on Eppler 387 drag.

Reference data: Williamson, McGranahan, Broughton, Deters, Brandt and Selig,
"Summary of Low-Speed Airfoil Data, Vol. 5" (UIUC, 2012), E387 drag polars
at Re = 1e5, 2e5, 3e5 and 4.6e5, and the UIUC E387 coordinates. They are
downloaded at run time, not stored here: UIUC's terms ask that their notices
accompany any redistribution.

Steps:
  python e387_calibrate.py fetch WORKDIR      # data + gapped coordinates
  # build e387_polar.cpp against the library, then for each Re and L:
  #   e387_polar WORKDIR/e387_gap.dat RE L 160 NCRIT > data/e387/PREFIX_RE_L.txt
  python e387_calibrate.py compare WORKDIR PREFIX L1 L2 ...
  python e387_calibrate.py docs WORKDIR ../../../docs/src/data/transition.json

The docs step writes the figures on the Transition validation page: the Cd
error at matched Cl (derived from, not a copy of, the UIUC data), the
summary table, and the dCd/dalpha sweeps in data/sweep_*_2412_1e6.txt.

The polars in data/e387 were made this way: prefix pol (ncrit 9) and pol8
(ncrit 8), 160 panels. Coordinates are reversed to lower-surface-first and
the trailing edge is opened to 0.25% chord (issue #18: closed trailing edges
fail). Above 160 panels the coarse E387 leading edge breaks the boundary-
layer initialization in Scirocco and mfoil alike.
"""
import os, re, sys, urllib.request
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
URLS = {'e387_drag.dat': 'https://m-selig.ae.illinois.edu/pd/pub/lsat/vol5/e387_drag.dat',
        'e387.dat': 'https://m-selig.ae.illinois.edu/ads/coord/e387.dat'}
TAGS = {1e5: '1e5', 2e5: '2e5', 3e5: '3e5', 4.6e5: '4.6e5'}


def fetch(work):
    os.makedirs(work, exist_ok=True)
    for name, url in URLS.items():
        urllib.request.urlretrieve(url, os.path.join(work, name))
    x = np.loadtxt(os.path.join(work, 'e387.dat'), skiprows=1)[::-1]  # lower first
    ile = np.argmin(x[:, 0]); g = 0.0025
    z = x[:, 1].copy()
    z[:ile] -= 0.5*g*x[:ile, 0]
    z[ile:] += 0.5*g*x[ile:, 0]
    np.savetxt(os.path.join(work, 'e387_gap.dat'), np.c_[x[:, 0], z],
               header='E387 TE gap 0.25%', comments='', fmt='%.6f')


def experiment(work):
    txt = open(os.path.join(work, 'e387_drag.dat')).read()
    blocks = {}
    pat = (r'Average Reynolds #:\s*([\d.]+)\s*Number of angles of attack:\s*(\d+)'
           r'\s*alpha[^\n]*\n((?:\s*-?[\d.]+(?:\s+-?[\d.]+)+\n){1,})')
    for m in re.finditer(pat, txt):
        rows = [[float(v) for v in l.split()[:3]] for l in m.group(3).strip().split('\n')]
        blocks[float(m.group(1))] = np.array(rows)
    return {Re: blocks[min(blocks, key=lambda r: abs(r - Re))] for Re in TAGS}


def errors(exp, prefix, Re, L):
    """Converged points with 0.1 <= Cl <= 1.0: (alpha, Cl, % Cd error)."""
    a, cl, cd = exp[Re].T
    sel = (cl >= 0.05) & (cl <= 1.05) & (a > -3.5)
    o = np.argsort(cl[sel])
    d = np.loadtxt(os.path.join(HERE, 'data', 'e387', '%s_%s_%d.txt' % (prefix, TAGS[Re], L)))
    d = d[(d[:, 3] == 1) & (d[:, 1] >= 0.1) & (d[:, 1] <= 1.0)]
    ref = np.interp(d[:, 1], cl[sel][o], cd[sel][o])
    return d[:, 0], d[:, 1], 100*(d[:, 2] - ref)/ref


def compare(work, prefix, Ls):
    """Cd error at matched Cl over 0.1 <= Cl <= 1.0 (attached flow)."""
    exp = experiment(work)
    print('Cd error vs UIUC E387: mean signed % / rms % / points')
    print('%-8s' % 'Re' + ''.join('%17s' % ('cut cell' if L == 0 else 'L=%d' % L) for L in Ls))
    allerr = {L: [] for L in Ls}
    for Re in TAGS:
        row = '%-8.1e' % Re
        for L in Ls:
            err = errors(exp, prefix, Re, L)[2]
            allerr[L] += list(err)
            row += '   %+5.1f /%5.1f /%2d' % (err.mean(), np.sqrt(np.mean(err**2)), len(err))
        print(row)
    print('%-8s' % 'all' + ''.join('   %+5.1f /%5.1f /%2d' % (
        np.mean(allerr[L]), np.sqrt(np.mean(np.square(allerr[L]))), len(allerr[L])) for L in Ls))


def docs(work, out):
    import json
    exp = experiment(work)
    r = lambda v, n: [round(float(x), n) for x in v]
    polars = []
    for Re in TAGS:
        ent = {'Re': Re}
        for key, L in (('cutcell', 0), ('captured', 60)):
            _, cl, err = errors(exp, 'pol', Re, L)
            ent[key] = {'cl': r(cl, 4), 'err': r(err, 2)}
        polars.append(ent)
    summary = []
    for prefix, ncrit in (('pol', 9), ('pol8', 8)):
        for L in (0, 40, 60, 80, 120, 240):
            if not os.path.exists(os.path.join(HERE, 'data', 'e387', '%s_1e5_%d.txt' % (prefix, L))):
                continue
            e = np.concatenate([errors(exp, prefix, Re, L)[2] for Re in TAGS])
            summary.append({'ncrit': ncrit, 'L': L, 'mean': round(float(e.mean()), 1),
                            'rms': round(float(np.sqrt(np.mean(e**2))), 1), 'n': int(e.size)})
    sweep = {}
    for key in ('cutcell', 'captured'):
        d = np.loadtxt(os.path.join(HERE, 'data', 'sweep_%s_2412_1e6.txt' % key))
        sweep[key] = {'alpha': r(d[:, 0], 3), 'dcd': r(d[:, 2], 8)}
    with open(out, 'w') as f:
        json.dump({'polars': polars, 'summary': summary, 'sweep': sweep}, f, separators=(',', ':'))


if __name__ == '__main__':
    if sys.argv[1] == 'fetch':
        fetch(sys.argv[2])
    elif sys.argv[1] == 'docs':
        docs(sys.argv[2], sys.argv[3])
    else:
        compare(sys.argv[2], sys.argv[3], [int(v) for v in sys.argv[4:]])
