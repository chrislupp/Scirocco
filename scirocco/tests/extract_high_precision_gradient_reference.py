#!/usr/bin/env python3
"""Full Python mfoil gradient reference in arbitrary precision.

Optional regeneration dependencies: numpy scipy matplotlib mpmath python-flint
(gmpy2 speeds up mpmath). Normal tests read the checked-in text fixture.
Run with MP_DPS=80 MP_STEPS=1e-14,1e-18 OPENBLAS_NUM_THREADS=1 for a
second precision/step check; the default is 60 digits and 1e-10,1e-14.

This evaluates Python's original panel closed forms, NOT the new quadrature
or C++ data. NumPy object arrays hold mpmath numbers; Flint solves the dense
linear systems at extra precision. The double solution supplies only a Newton
initial guess. Tighten the wake-spacing and transition subsolves, then solve
the complete coupled equations with residual below 10**(-MP_DPS+15).
The ten-pass not-a-knot arc-length spline is reproduced below because SciPy
cannot interpolate arbitrary-precision arrays. Keep the original model's
binary floating-point constants and quadrature data. Pressure outputs are
normalized once, as in extract_remaining_gradient_reference.py.

The central-difference step/precision study measures this smooth solution
branch. It is not interval certification of the entire nonlinear algorithm.
"""
import copy
import hashlib
import os
import time
import types
from pathlib import Path

import flint
import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
mp.mp.dps = int(os.environ.get("MP_DPS", "60"))
if mp.mp.dps < 50:
    raise ValueError("Use at least 50 decimal digits for this reference")
flint.ctx.dps = mp.mp.dps + 15
start = time.monotonic()


def note(*a):
    print(round(time.monotonic() - start, 2), *a, flush=True)


def scalar(x):
    if isinstance(x, np.ndarray) and x.ndim == 0:
        x = x.item()
    if isinstance(x, (np.integer, int)):
        return mp.mpf(int(x))
    return mp.mpf(x)


def arr(x):
    return np.vectorize(scalar, otypes=[object])(np.asarray(x))


def linear_solve(A, b):
    A = np.asarray(A)
    b = np.asarray(b)
    isvec = b.ndim == 1
    if isvec:
        b = b[:, None]

    def toarb(x):
        return flint.arb(mp.nstr(scalar(x), mp.mp.dps + 3))

    aa = flint.arb_mat([[toarb(x) for x in row] for row in A])
    bb = flint.arb_mat([[toarb(x) for x in row] for row in b])
    z = aa.solve(bb)
    out = np.array(
        [
            [
                mp.mpf(z[i, j].mid().str(mp.mp.dps + 3, radius=False))
                for j in range(b.shape[1])
            ]
            for i in range(A.shape[0])
        ],
        dtype=object,
    )
    return out[:, 0] if isvec else out


class SparseMP:
    __array_priority__ = 1000

    def __init__(self, value):
        self.a = (
            np.full(value, mp.mpf(0), dtype=object)
            if isinstance(value, tuple)
            else arr(value)
        )

    @property
    def shape(self):
        return self.a.shape

    def __getitem__(self, k):
        v = self.a[k]
        if isinstance(k, tuple) and len(k) == 2 and np.ndim(v) == 1:
            return v[:, None] if np.isscalar(k[1]) else v[None, :]
        return v

    def __setitem__(self, k, v):
        if isinstance(v, SparseMP):
            v = v.a
        if self.a[k].ndim == 1:
            v = np.asarray(v).ravel()
        self.a[k] = v

    def __matmul__(self, b):
        return self.a @ (b.a if isinstance(b, SparseMP) else b)

    def __rmatmul__(self, b):
        return b @ self.a

    def __imul__(self, v):
        self.a *= v
        return self

    def tocsr(self):
        return self

    def toarray(self):
        return self.a.copy()


class NumpyMP:
    pi = mp.pi
    linalg = types.SimpleNamespace(
        solve=linear_solve,
        norm=lambda x, *args: mp.sqrt(sum(v * v for v in np.asarray(x).ravel())),
    )

    def __getattr__(self, k):
        return getattr(np, k)

    def zeros(self, shape, dtype=None):
        return (
            np.zeros(shape, dtype=dtype)
            if dtype in (int, bool)
            else np.full(shape, mp.mpf(0), dtype=object)
        )

    def ones(self, shape, dtype=None):
        return (
            np.ones(shape, dtype=dtype)
            if dtype in (int, bool)
            else np.full(shape, mp.mpf(1), dtype=object)
        )

    def array(self, x, dtype=None, **kw):
        return np.array(x, dtype=dtype, **kw) if dtype in (int, bool) else arr(x)

    def eye(self, *a, **kw):
        return arr(np.eye(*a, **kw))


for n, fn in [
    ("sqrt", mp.sqrt),
    ("log", mp.log),
    ("log10", mp.log10),
    ("exp", mp.exp),
    ("tanh", mp.tanh),
    ("sin", mp.sin),
    ("cos", mp.cos),
    ("ceil", mp.ceil),
    ("arctan2", mp.atan2),
]:

    def mapped(*args, _fn=fn):
        if all(np.ndim(a) == 0 for a in args):
            return _fn(*(scalar(a) for a in args))
        return np.vectorize(_fn, otypes=[object])(*args)

    setattr(NumpyMP, n, staticmethod(mapped))
source = (ROOT / "mfoil/mfoil.py").read_text()
# Solve the same implicit equations more tightly for this precision study.
for old in ("if (abs(Rxt) < 1e-14): break", "if (abs(dr)<1e-14): break"):
    assert source.count(old) == 1, "Python reference stopping rule changed"
    source = source.replace(old, old[: old.index("1e-14")] + "_reference_tol): break")
m = types.ModuleType("mpfoil")
m.__dict__["_reference_tol"] = mp.mpf(10) ** (-mp.mp.dps + 10)
exec(compile(source, str(ROOT / "mfoil/mfoil.py"), "exec"), m.__dict__)
base = m.mfoil(naca="0012", npanel=100)
base.param.verb = 0
base.param.doplot = False
x = 2 * base.foil.x.copy()
base.foil.x = x
base.geom.chord = 2.0
base.foil.s = m.spline2d(x)["X"].x.copy()
base.setoper(alpha=3, Re=1e6, Ma=0.3)
base.solve()
assert base.glob.conv
note("double initialized")


# Convert the solver inputs and state, leaving flags/index arrays integral.
def convert(v):
    if isinstance(v, np.ndarray):
        return v.copy() if v.dtype.kind in "biu" else arr(v)
    if isinstance(v, (float, np.floating)):
        return scalar(v)
    if isinstance(v, list):
        return [convert(a) for a in v]
    if isinstance(v, tuple):
        return tuple(convert(a) for a in v)
    if hasattr(v, "toarray"):
        return SparseMP(v.toarray())
    if hasattr(v, "__dict__"):
        z = copy.copy(v)
        for k, a in vars(v).items():
            setattr(z, k, convert(a))
        return z
    return copy.deepcopy(v)


base = convert(base)
m.np = NumpyMP()
m.sparse = types.SimpleNamespace(
    lil_matrix=SparseMP, spdiags=lambda data, *args: SparseMP(np.diag(data))
)
# Use the original closed forms at high precision, independent of the newly
# introduced double-precision quadrature evaluators.
m.panel_linvortex_stream_far = lambda *args: None
m.panel_source_stream_far = lambda *args: None
m.vprint = lambda *args: None
m._reference_tol = mp.mpf(10) ** (-mp.mp.dps + 10)
q, w = map(arr, m.quadseg())


def geometry(M):
    X = M.foil.x
    n = X.shape[1]
    s = [mp.mpf(0)]
    for i in range(1, n):
        s.append(s[-1] + mp.sqrt(sum((X[d, i] - X[d, i - 1]) ** 2 for d in (0, 1))))

    def spline(y):
        h = [s[i + 1] - s[i] for i in range(n - 1)]
        sec = [(y[i + 1] - y[i]) / h[i] for i in range(n - 1)]
        lower = [mp.mpf(0)] * n
        diag = lower.copy()
        upper = lower.copy()
        b = lower.copy()
        for i in range(1, n - 1):
            lower[i] = h[i]
            diag[i] = 2 * (h[i - 1] + h[i])
            upper[i] = h[i - 1]
            b[i] = 3 * (h[i] * sec[i - 1] + h[i - 1] * sec[i])
        d = h[0] + h[1]
        diag[0] = h[1]
        upper[0] = d
        b[0] = ((h[0] + 2 * d) * h[1] * sec[0] + h[0] ** 2 * sec[1]) / d
        d = h[-2] + h[-1]
        lower[-1] = d
        diag[-1] = h[-2]
        b[-1] = (h[-1] ** 2 * sec[-2] + (2 * d + h[-1]) * h[-2] * sec[-1]) / d
        for i in range(1, n):
            f = lower[i] / diag[i - 1]
            diag[i] -= f * upper[i - 1]
            b[i] -= f * b[i - 1]
        slopes = [mp.mpf(0)] * n
        slopes[-1] = b[-1] / diag[-1]
        for i in range(n - 2, -1, -1):
            slopes[i] = (b[i] - upper[i] * slopes[i + 1]) / diag[i]
        out = []
        for i in range(n - 1):
            t = (slopes[i] + slopes[i + 1] - 2 * sec[i]) / h[i]
            out.append([t / h[i], (sec[i] - slopes[i]) / h[i] - t, slopes[i]])
        return out

    for _ in range(10):
        pp = [spline(y) for y in X]
        new = [mp.mpf(0)]
        for i in range(n - 1):
            h = s[i + 1] - s[i]
            integral = 0
            for qq, ww in zip(q, w):
                t = qq * h
                v = [(3 * c[i][0] * t + 2 * c[i][1]) * t + c[i][2] for c in pp]
                integral += ww * mp.sqrt(v[0] ** 2 + v[1] ** 2)
            new.append(new[-1] + h * integral)
        s = new
    M.foil.s = np.array(s, dtype=object)


def residual(M):
    m.stagpoint_move(M)
    m.set_wake_gap(M)
    m.build_glob_sys(M)
    ue = M.glob.U[3]
    ds = M.glob.U[1]
    ueinv = m.get_ueinv(M)
    n = len(ue)
    R = np.r_[M.glob.R, ue - ueinv - M.vsol.ue_m @ (ds * ue)]
    A = m.np.zeros((4 * n, 4 * n))
    A[: 3 * n] = M.glob.R_U.a
    A[3 * n :, 3::4] = m.np.eye(n) - M.vsol.ue_m * ds[None, :]
    A[3 * n :, 1::4] = -M.vsol.ue_m * ue[None, :]
    return R, A


def run(delta):
    M = copy.deepcopy(base)
    M.foil.x[0, 50] += delta
    geometry(M)
    note("geometry", str(delta))
    m.init_thermo(M)
    m.build_gamma(M, M.oper.alpha)
    note("gamma")
    m.build_wake(M)
    note("wake")
    m.stagpoint_move(M)
    m.set_wake_gap(M)
    m.calc_ue_m(M)
    note("coupling")
    for it in range(8):
        R, A = residual(M)
        rn = max(abs(r) for r in R)
        note("Newton", it, mp.nstr(rn, 8))
        if rn < mp.mpf(10) ** (-mp.mp.dps + 15):
            break
        du = linear_solve(A, -R).reshape(4, -1, order="F")
        M.glob.U += du
    else:
        raise RuntimeError("MP Newton did not converge")
    chord = M.geom.chord
    M.geom.chord = mp.mpf(1)
    m.calc_force(M)
    M.geom.chord = chord
    values = [M.post.cl / chord, M.post.cd, M.post.cm / chord**2, M.post.cdpi / chord]
    note("outputs", *[mp.nstr(v, 30) for v in values])
    return values


if __name__ == "__main__":
    results = []
    for h in [mp.mpf(s) for s in os.environ.get("MP_STEPS", "1e-10,1e-14").split(",")]:
        plus = run(h)
        minus = run(-h)
        deriv = [(a - b) / (2 * h) for a, b in zip(plus, minus)]
        results.append((h, deriv))
        note("central", str(h), *[mp.nstr(v, 40) for v in deriv])
    if len(results) < 2:
        raise RuntimeError("Use at least two steps to check convergence")
    assert max(abs(a - b) for a, b in zip(results[-1][1], results[-2][1])) < mp.mpf(
        "2e-16"
    )
    output = Path(
        os.environ.get(
            "MP_OUTPUT",
            Path(__file__).with_name("high_precision_gradient_reference.txt"),
        )
    )
    with output.open("w") as f:
        f.write(
            "# NACA 0012; 100 panels; chord 2; alpha 3; Re 1e6; Ma .3; spline arc labels\n"
        )
        f.write(
            "# d(Cl,Cd,Cm,Cdpi)/d(x at node 50), original Python panel closed forms\n"
        )
        f.write(
            "# mfoil.py SHA256 "
            + hashlib.sha256((ROOT / "mfoil/mfoil.py").read_bytes()).hexdigest()
            + "\n"
        )
        f.write(
            "# decimal digits "
            + str(mp.mp.dps)
            + "; central steps "
            + ",".join(str(h) for h, _ in results)
            + "\n"
        )
        f.write(" ".join(mp.nstr(v, 40) for v in results[-1][1]) + "\n")
    note("wrote", output)
