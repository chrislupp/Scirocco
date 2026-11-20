#!/usr/bin/env python3
"""Data for the documentation's verification and validation pages.

    vv.py codes       DRIVER OUT.json   Scirocco, XFOIL and mfoil
    vv.py experiments DRIVER OUT.json   Scirocco and XFOIL against wind-tunnel data

DRIVER is vv_polar built against the library (see README.md). Every code runs
on the same panel nodes: mfoil panels the airfoil, Scirocco takes the nodes,
arc lengths and tangents unchanged, and XFOIL loads the nodes without
repaneling (PCOP). Scirocco and XFOIL sweep each polar in order, as a user
would, so every point warm-starts from the previous one; mfoil solves each
point from a cold start.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.io import FortranFile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "mfoil"))
import mfoil  # noqa: E402

XFOIL = os.environ.get("XFOIL", "xfoil")
NPANEL = 199
NCRIT = 9.0


def arange(a0, a1, step):
    return [round(a, 4) for a in np.arange(a0, a1 + 1e-9, step)]


# Code-to-code cases: (name, naca, Re, Mach, alphas, distribution alphas).
CODE_CASES = [
    ("inviscid_2412", "2412", 0, 0.0, arange(-4, 12, 2), [4.0]),
    ("0012_re1e6", "0012", 1e6, 0.0, arange(0, 15, 1), [8.0]),
    ("2412_re3e6", "2412", 3e6, 0.0, arange(-4, 14, 1), [4.0]),
    ("4412_re2e5", "4412", 2e5, 0.0, arange(-2, 12, 1), [4.0]),
    ("2412_re6e6_m04", "2412", 6e6, 0.4, arange(-2, 6, 1), [2.0]),
]

# Validation cases against published measurements.
EXPERIMENT_CASES = [
    ("0012_re6e6", "0012", 6e6, 0.15, arange(0, 19, 1), []),
    ("0012_re288e4", "0012", 2.88e6, 0.0, arange(0, 15, 1), [0.0, 10.0, 15.0]),
    ("4412_re152e4", "4412", 1.52e6, 0.08, arange(0, 13, 1) + [13.87], [13.87]),
]


def sweeps(alphas):
    """Split a polar into two sweeps that both start at the angle nearest 0."""
    alphas = sorted(alphas)
    k = int(np.argmin(np.abs(alphas)))
    up, down = alphas[k:], alphas[k::-1]
    return [s for s in (up, down) if len(s) > 1 or s is up]


def panel(naca, workdir):
    m = mfoil.mfoil(naca=naca, npanel=NPANEL)
    nodes = np.column_stack((m.foil.x.T, m.foil.s, m.foil.t.T))
    np.savetxt(workdir / "nodes.txt", nodes, fmt="%.17g")
    # XFOIL orders the nodes from the upper trailing edge.
    np.savetxt(workdir / "airfoil.dat", m.foil.x.T[::-1], fmt="%.17g",
               header=f"NACA{naca}", comments="")
    return m.foil.x


def run_scirocco(driver, workdir, re, ma, alphas, dist):
    points, dists = {}, {}
    for k, sweep in enumerate(sweeps(alphas)):
        out = workdir / f"scirocco_{k}.json"
        subprocess.run(
            [str(driver), str(workdir / "nodes.txt"), str(re), str(ma), str(NCRIT),
             ",".join(map(str, sweep)),
             ",".join(str(a) for a in dist if a in sweep) or "-", str(out)],
            check=True, timeout=600)
        data = json.loads(out.read_text())
        for p in data["points"]:
            points.setdefault(p["alpha"], p)
        for d in data["distributions"]:
            dists.setdefault(d["alpha"], d)
    return points, dists


def run_xfoil(workdir, re, ma, alphas, dist, n):
    points, dists = {}, {}
    for k, sweep in enumerate(sweeps(alphas)):
        tag = f"xfoil_{k}"
        dump = workdir / f"{tag}.bin"
        dump.unlink(missing_ok=True)
        cmds = ["PLOP", "G F", "", "LOAD airfoil.dat", "PCOP", "OPER",
                f"MACH {ma}", "ITER 100", "VPAR", f"N {NCRIT:g}", "XTR 1 1", ""]
        if re:
            cmds.append(f"VISC {re:g}")
        cmds += ["PACC", "", dump.name]
        cmds += [f"ALFA {a}" for a in sweep]
        cmds += ["", "QUIT", ""]
        subprocess.run([XFOIL], input="\n".join(cmds), text=True, cwd=workdir,
                       capture_output=True, timeout=600)
        with FortranFile(dump, "r") as f:
            f.read_record(np.uint8)
            f.read_reals(np.float64)
            f.read_ints(np.int32)
            info = f.read_ints(np.int32)
            assert info[3] == n, "XFOIL repaneled the airfoil"
            nodes = f.read_reals(np.float64).reshape(-1, 2)
            while True:
                try:
                    v = f.read_reals(np.float64)
                except Exception:
                    break
                alpha = round(float(v[0]), 4)
                p = dict(alpha=alpha, converged=True, cl=v[1], cd=v[2], cm=v[4])
                if re:
                    p.update(xtr_upper=v[5], xtr_lower=v[6])
                bl = f.read_ints(np.int32)
                upper = f.read_reals(np.float64).reshape(-1, 6)
                lower = f.read_reals(np.float64).reshape(-1, 6)
                points.setdefault(alpha, p)
                if alpha in dist and alpha not in dists:
                    # Index 0 of each side is a synthetic stagnation state;
                    # wake points follow the surface nodes.
                    s = np.vstack((upper[1:bl[2]][::-1], lower[1:bl[3]]))[::-1]
                    d = dict(alpha=alpha, x=nodes[::-1, 0], y=nodes[::-1, 1],
                             cp=s[:, 1])
                    if re:
                        d.update(cf=s[:, 4], theta=s[:, 2], dstar=s[:, 3])
                    dists[alpha] = d
    for a in alphas:
        points.setdefault(a, dict(alpha=a, converged=False))
    return points, dists


def mfoil_point(args):
    naca, re, ma, alpha, want_dist = args
    m = mfoil.mfoil(naca=naca, npanel=NPANEL)
    m.param.verb = 0
    m.param.doplot = False
    m.param.ncrit = NCRIT
    m.setoper(alpha=alpha, Ma=ma)
    if re:
        m.setoper(Re=re)
    try:
        m.solve()
    except Exception:
        return dict(alpha=alpha, converged=False), None
    ok = bool(m.glob.conv) if re else True
    p = dict(alpha=alpha, converged=ok)
    if not ok:
        return p, None
    p.update(cl=m.post.cl, cd=m.post.cd, cm=m.post.cm)
    if re:
        # Xt holds (xi, x) per surface; x stays 0 if the surface never
        # transitions, which XFOIL reports as 1.
        p.update(xtr_upper=m.vsol.Xt[1, 1] or 1.0,
                 xtr_lower=m.vsol.Xt[0, 1] or 1.0)
    d = None
    if want_dist:
        n = m.foil.N
        d = dict(alpha=alpha, x=m.foil.x[0], y=m.foil.x[1], cp=m.post.cp[:n])
        if re:
            mfoil.get_distributions(m)
            d.update(cf=m.post.cf[:n], theta=m.post.th[:n], dstar=m.post.ds[:n])
    return p, d


def run_mfoil(pool, naca, re, ma, alphas, dist):
    results = list(pool.map(mfoil_point,
                            [(naca, re, ma, a, a in dist) for a in alphas]))
    points = {p["alpha"]: p for p, _ in results}
    dists = {d["alpha"]: d for _, d in results if d}
    return points, dists


def clean(obj, digits=6):
    """Round floats and convert NumPy values for a compact JSON file."""
    if isinstance(obj, dict):
        return {k: clean(v, digits) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, np.ndarray)):
        return [clean(v, digits) for v in obj]
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        v = float(obj)
        return None if not np.isfinite(v) else float(f"{v:.{digits}g}")
    return obj


def run_case(case, driver, codes, pool):
    name, naca, re, ma, alphas, dist = case
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        x = panel(naca, workdir)
        out = dict(name=name, airfoil=f"NACA {naca}", re=re, mach=ma,
                   ncrit=NCRIT, panels=NPANEL, alphas=alphas,
                   coords=dict(x=x[0], y=x[1]))
        runs = {}
        if "scirocco" in codes:
            runs["scirocco"] = run_scirocco(driver, workdir, re, ma, alphas, dist)
        if "xfoil" in codes:
            runs["xfoil"] = run_xfoil(workdir, re, ma, alphas, dist, x.shape[1])
        if "mfoil" in codes:
            runs["mfoil"] = run_mfoil(pool, naca, re, ma, alphas, dist)
    for code, (points, dists) in runs.items():
        out[code] = dict(points=[points[a] for a in alphas],
                         distributions=[dists[a] for a in dist if a in dists])
        bad = [a for a in alphas if not points[a]["converged"]]
        print(f"{name:16s} {code:9s} unconverged: {bad}", flush=True)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("mode", choices=("codes", "experiments"))
    parser.add_argument("driver", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    args = parser.parse_args()
    driver = args.driver.resolve()
    with ProcessPoolExecutor(args.workers) as pool:
        if args.mode == "codes":
            cases = [run_case(c, driver, ("scirocco", "xfoil", "mfoil"), pool)
                     for c in CODE_CASES]
            data = dict(cases=cases)
        else:
            cases = [run_case(c, driver, ("scirocco", "xfoil"), pool)
                     for c in EXPERIMENT_CASES]
            data = dict(cases=cases, measured=load_measurements())
    args.output.write_text(json.dumps(clean(data), separators=(",", ":")) + "\n")


def load_measurements():
    """Tabulated measurements from experiments/ (see its README)."""
    def table(name):
        rows, zone = {}, None
        for line in (HERE / "experiments" / name).read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("variables"):
                continue
            if line.lower().startswith("zone"):
                zone = line.split('"')[1]
                continue
            rows.setdefault(zone, []).append([float(v) for v in line.split()])
        return {k: np.array(v) for k, v in rows.items()}

    def only(name):
        (zone,) = table(name).values()
        return zone

    abbott_cl = only("0012.abbottdata.cl.dat")
    abbott_cd = only("0012.abbottdata.cd.dat")
    ladson = table("CLCD_Ladson_expdata.dat")
    gregory_cl = only("CL_Gregory_expdata.dat")
    gregory_cp = table("CP_Gregory_expdata.dat")
    coles = only("naca4412.cp.expt.dat")
    return dict(
        abbott_cl=dict(alpha=abbott_cl[:, 0], cl=abbott_cl[:, 1]),
        abbott_cd=dict(cl=abbott_cd[:, 0], cd=abbott_cd[:, 1]),
        ladson={z: dict(alpha=v[:, 0], cl=v[:, 1], cd=v[:, 2])
                for z, v in ladson.items()},
        gregory_cl=dict(alpha=gregory_cl[:, 0], cl=gregory_cl[:, 1]),
        gregory_cp={z: dict(x=v[:, 0], cp=v[:, 1]) for z, v in gregory_cp.items()},
        coles_wadcock_cp=dict(x=coles[:, 0], y=coles[:, 1], cp=coles[:, 2]),
    )


if __name__ == "__main__":
    main()
