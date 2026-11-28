#!/usr/bin/env python3
"""Data for the documentation's verification and validation pages.

    vv.py codes       DRIVER OUT.json   Scirocco (mfoil's closure), XFOIL, mfoil
                                        and flexfoil
    vv.py experiments DRIVER OUT.json   Scirocco and XFOIL against wind-tunnel data
    vv.py tripped     DRIVER OUT.json   Scirocco (both transition models, and
                                        without the skin-friction floor) and
                                        XFOIL against tripped low-Reynolds-
                                        number polars (UIUC)

DRIVER is vv_polar built against the library (see README.md). Every code runs
on the same panel nodes: mfoil panels the airfoil, Scirocco takes the nodes,
arc lengths and tangents unchanged, and XFOIL loads the nodes without
repaneling (PCOP), as does flexfoil. Scirocco and XFOIL sweep each polar in
order, as a user would, so every point warm-starts from the previous one;
mfoil and flexfoil solve each point from a cold start.

tripped downloads the UIUC measurements and coordinates into --cache (they
are GPL-licensed and not kept in this repository). XFOIL panels each airfoil
(160 panels) and both codes use those nodes.
"""
import argparse
import json
import os
import subprocess
import sys
import re as regex
import tempfile
import urllib.request
import zipfile
from concurrent.futures import ProcessPoolExecutor
from importlib import metadata
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

# Validation cases against published measurements. A seventh entry is forced
# transition, (x/c upper, x/c lower); Ladson's roughness strips were at 5%
# chord on both surfaces.
EXPERIMENT_CASES = [
    ("0012_re6e6", "0012", 6e6, 0.15, arange(0, 19, 1), []),
    ("0012_re6e6_trip", "0012", 6e6, 0.15, arange(0, 19, 1), [], (0.05, 0.05)),
    ("0012_re288e4", "0012", 2.88e6, 0.0, arange(0, 15, 1), [0.0, 10.0, 15.0]),
    ("4412_re152e4", "4412", 1.52e6, 0.08, arange(0, 13, 1) + [13.87], [13.87]),
]


# Tripped polars from Selig & McGranahan, "Wind Tunnel Aerodynamic Tests of
# Six Airfoils for Use on Small Wind Turbines", NREL/SR-500-34515, 2004 (UIUC
# Low-Speed Airfoil Tests, vol. 4): zigzag trips at 2% chord on the upper
# and 5% on the lower surface. (name, label, coordinates, drag file); the
# other four airfoils of the report are not run.
UIUC = "https://m-selig.ae.illinois.edu/"
UIUC_DATA = "pd/pub/lsat/volume04.zip"
UIUC_TRIP = (0.02, 0.05)
UIUC_PANELS = 160
UIUC_TE_GAP = 0.0025  # the solver needs an open trailing edge
TRIPPED_CASES = [
    ("e387", "Eppler 387", "ads/coord/e387.dat", "volume04/e387_tf_drg.txt"),
    ("sd2030", "SD2030", "ads/coord/sd2030.dat", "volume04/sd2030_tf_drg.txt"),
]
TRIPPED_ALPHAS = arange(-3, 9, 0.5)
# (key, vv_polar flags)
TRIPPED_RUNS = [("scirocco", []), ("scirocco_captured", ["captured"]),
                ("scirocco_nofloor", ["mfoil"])]


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


def run_scirocco(driver, workdir, re, ma, alphas, dist, trip, flags=()):
    """flags: vv_polar's trailing arguments ("mfoil", "captured")."""
    points, dists = {}, {}
    for k, sweep in enumerate(sweeps(alphas)):
        out = workdir / f"scirocco_{'_'.join(flags)}_{k}.json"
        subprocess.run(
            [str(driver), str(workdir / "nodes.txt"), str(re), str(ma), str(NCRIT),
             ",".join(map(str, sweep)),
             ",".join(str(a) for a in dist if a in sweep) or "-", str(out)]
            + [f"{trip[0]},{trip[1]}" if trip else "-"] + list(flags),
            check=True, timeout=600)
        data = json.loads(out.read_text())
        for p in data["points"]:
            points.setdefault(p["alpha"], p)
        for d in data["distributions"]:
            dists.setdefault(d["alpha"], d)
    return points, dists


def run_xfoil(workdir, re, ma, alphas, dist, n, trip):
    xtr = "XTR {:g} {:g}".format(*trip) if trip else "XTR 1 1"
    points, dists = {}, {}
    for k, sweep in enumerate(sweeps(alphas)):
        tag = f"xfoil_{k}"
        dump = workdir / f"{tag}.bin"
        dump.unlink(missing_ok=True)
        cmds = ["PLOP", "G F", "", "LOAD airfoil.dat", "PCOP", "OPER",
                f"MACH {ma}", "ITER 100", "VPAR", f"N {NCRIT:g}", xtr, ""]
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


def flexfoil_coords(x):
    """flexfoil's flat coordinate list, from the upper trailing edge, as XFOIL."""
    return [float(v) for v in x.T[::-1].ravel()]


def flexfoil_dist(b, alpha, ma):
    """Surface distribution from flexfoil's boundary-layer stations.

    Station 0 of each side is the stagnation point; the rest are the panel
    nodes. Cp follows from the edge velocity with XFOIL's Karman-Tsien
    correction (CPCALC). flexfoil normalizes Cf by the local edge velocity,
    so Cf * ue^2 is XFOIL's Cf, based on the freestream; that holds only for
    incompressible flow, so Cf is left out when M > 0.
    """
    def side(key):
        return np.concatenate((np.array(b[f"{key}_lower"][1:])[::-1],
                               b[f"{key}_upper"][1:]))
    ue = side("ue")
    cpinc = 1.0 - ue**2
    beta = np.sqrt(1.0 - ma**2)
    cp = cpinc / (beta + 0.5 * ma**2 / (1.0 + beta) * cpinc)
    d = dict(alpha=alpha, x=side("x"), cp=cp, theta=side("theta"),
             dstar=side("delta_star"))
    if ma == 0:
        d["cf"] = side("cf") * ue**2
    return d


def flexfoil_point(args):
    flat, re, ma, alpha, want_dist = args
    from flexfoil import _rustfoil as ff
    failed = dict(alpha=alpha, converged=False), None
    if not re:
        r = ff.analyze_inviscid(flat, alpha)
        if not r["success"]:
            return failed
        p = dict(alpha=alpha, converged=True, cl=r["cl"], cd=0.0, cm=r["cm"])
        # The inviscid solver reports Cp at panel midpoints.
        d = dict(alpha=alpha, x=r["cp_x"], cp=r["cp"]) if want_dist else None
        return p, d
    r = ff.analyze_faithful(flat, alpha, reynolds=re, mach=ma, ncrit=NCRIT,
                            max_iterations=100)
    if not (r["success"] and r["converged"]):
        return failed
    p = dict(alpha=alpha, converged=True, cl=r["cl"], cd=r["cd"], cm=r["cm"],
             xtr_upper=r["x_tr_upper"], xtr_lower=r["x_tr_lower"])
    d = None
    if want_dist:
        # analyze_faithful returns no distributions; solve the point again.
        b = ff.get_bl_distribution(flat, alpha, reynolds=re, mach=ma,
                                   ncrit=NCRIT, max_iterations=100)
        if b["success"] and b["converged"]:
            d = flexfoil_dist(b, alpha, ma)
    return p, d


def run_flexfoil(pool, x, re, ma, alphas, dist):
    flat = flexfoil_coords(x)
    results = list(pool.map(flexfoil_point,
                            [(flat, re, ma, a, a in dist) for a in alphas]))
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
    name, naca, re, ma, alphas, dist, *rest = case
    trip = rest[0] if rest else None
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        x = panel(naca, workdir)
        out = dict(name=name, airfoil=f"NACA {naca}", re=re, mach=ma,
                   ncrit=NCRIT, panels=NPANEL, alphas=alphas,
                   coords=dict(x=x[0], y=x[1]))
        if trip:
            out.update(xtrip_upper=trip[0], xtrip_lower=trip[1])
        runs = {}
        if "scirocco" in codes:
            # Against mfoil, Scirocco runs mfoil's closure: without the
            # laminar floor on the turbulent skin friction that is its default.
            runs["scirocco"] = run_scirocco(driver, workdir, re, ma, alphas, dist,
                                            trip, ["mfoil"] * ("mfoil" in codes))
        if "scirocco_captured" in codes:
            runs["scirocco_captured"] = run_scirocco(
                driver, workdir, re, ma, alphas, dist, trip, ["captured"])
        if "xfoil" in codes:
            runs["xfoil"] = run_xfoil(workdir, re, ma, alphas, dist, x.shape[1],
                                      trip)
        if "mfoil" in codes:
            if trip:
                raise ValueError("mfoil has no forced transition")
            runs["mfoil"] = run_mfoil(pool, naca, re, ma, alphas, dist)
        if "flexfoil" in codes:
            runs["flexfoil"] = run_flexfoil(pool, x, re, ma, alphas, dist)
    for code, (points, dists) in runs.items():
        out[code] = dict(points=[points[a] for a in alphas],
                         distributions=[dists[a] for a in dist if a in dists])
        bad = [a for a in alphas if not points[a]["converged"]]
        print(f"{name:16s} {code:9s} unconverged: {bad}", flush=True)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("mode", choices=("codes", "experiments", "tripped"))
    parser.add_argument("driver", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--cache", type=Path, default=Path("uiuc_cache"),
                        help="tripped: where the UIUC downloads are kept")
    args = parser.parse_args()
    driver = args.driver.resolve()
    if args.mode == "tripped":
        data = dict(cases=run_tripped(driver, args.cache),
                    trip=dict(upper=UIUC_TRIP[0], lower=UIUC_TRIP[1]),
                    panels=UIUC_PANELS, ncrit=NCRIT, cl_range=[0.1, 1.0])
        args.output.write_text(
            json.dumps(clean(data), separators=(",", ":")) + "\n")
        return
    with ProcessPoolExecutor(args.workers) as pool:
        if args.mode == "codes":
            codes = ("scirocco", "xfoil", "mfoil", "flexfoil")
            cases = [run_case(c, driver, codes, pool) for c in CODE_CASES]
            data = dict(cases=cases,
                        versions=dict(flexfoil=metadata.version("flexfoil")))
        else:
            cases = [run_case(c, driver,
                              ("scirocco", "xfoil")
                              + ("scirocco_captured",) * (len(c) > 6), pool)
                     for c in EXPERIMENT_CASES]
            data = dict(cases=cases, measured=load_measurements())
    args.output.write_text(json.dumps(clean(data), separators=(",", ":")) + "\n")


def uiuc_fetch(cache, path):
    """Local copy of a UIUC file, downloaded on first use."""
    local = cache / Path(path).name
    if not local.exists():
        cache.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(UIUC + path, local)
    return local


def uiuc_polars(text):
    """{average Reynolds number: rows of alpha, Cl, Cd} of a UIUC drag file."""
    pattern = (r"Average Reynolds #:\s*([\d.]+)\s*Number of angles of attack:"
               r"\s*\d+\s*alpha[^\n]*\n((?:\s*-?[\d.]+(?:\s+-?[\d.]+)+\n)+)")
    return {float(m.group(1)):
            np.array([[float(v) for v in row.split()[:3]]
                      for row in m.group(2).strip().split("\n")])
            for m in regex.finditer(pattern, text)}


def uiuc_panel(coords, workdir):
    """Open the trailing edge and let XFOIL panel the airfoil.

    Writes airfoil.dat (XFOIL's order) and nodes.txt (x y, from the lower
    trailing edge, for vv_polar). The published coordinates are too coarse
    at the leading edge to panel directly.
    """
    x = np.loadtxt(coords, skiprows=1)  # from the upper trailing edge
    le = int(np.argmin(x[:, 0]))
    x[:le, 1] += 0.5 * UIUC_TE_GAP * x[:le, 0]
    x[le:, 1] -= 0.5 * UIUC_TE_GAP * x[le:, 0]
    np.savetxt(workdir / "raw.dat", x, fmt="%.7f", header="raw", comments="")
    cmds = ["PLOP", "G F", "", "LOAD raw.dat", "PPAR", f"N {UIUC_PANELS + 1}",
            "", "", "PSAV airfoil.dat", "QUIT", ""]
    subprocess.run([XFOIL], input="\n".join(cmds), text=True, cwd=workdir,
                   capture_output=True, timeout=60)
    nodes = np.loadtxt(workdir / "airfoil.dat")
    assert len(nodes) == UIUC_PANELS + 1, "XFOIL did not panel the airfoil"
    np.savetxt(workdir / "airfoil.dat", nodes, fmt="%.17g", header="airfoil",
               comments="")
    np.savetxt(workdir / "nodes.txt", nodes[::-1], fmt="%.17g")
    return nodes[::-1].T


def cd_error(measured, points, cl_lo=0.1, cl_hi=1.0):
    """Cd error [%] at the measured Cl-Cd curve, for converged points with
    cl_lo <= Cl <= cl_hi: mean, rms and number of points."""
    alpha, cl, cd = measured.T
    keep = (cl >= cl_lo - 0.05) & (cl <= cl_hi + 0.05) & (alpha > -3.5)
    order = np.argsort(cl[keep])
    pts = [p for p in points
           if p["converged"] and cl_lo <= p["cl"] <= cl_hi]
    if not pts:
        return dict(mean=None, rms=None, points=0)
    ref = np.interp([p["cl"] for p in pts], cl[keep][order], cd[keep][order])
    err = 100 * (np.array([p["cd"] for p in pts]) - ref) / ref
    return dict(mean=err.mean(), rms=np.sqrt((err**2).mean()), points=len(pts))


def run_tripped(driver, cache):
    archive = zipfile.ZipFile(uiuc_fetch(cache, UIUC_DATA))
    cases = []
    for name, label, coords, drag in TRIPPED_CASES:
        polars = uiuc_polars(
            archive.read(drag).decode("latin-1").replace("\r\n", "\n"))
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            x = uiuc_panel(uiuc_fetch(cache, coords), workdir)
            for re_measured, measured in sorted(polars.items()):
                re = round(re_measured, -4)
                out = dict(name=f"{name}_re{re:.0f}", airfoil=label, re=re,
                           re_measured=re_measured, alphas=TRIPPED_ALPHAS,
                           coords=dict(x=x[0], y=x[1]),
                           measured=dict(alpha=measured[:, 0],
                                         cl=measured[:, 1], cd=measured[:, 2]),
                           summary={})
                runs = {key: run_scirocco(driver, workdir, re, 0.0,
                                          TRIPPED_ALPHAS, [], UIUC_TRIP, flags)
                        for key, flags in TRIPPED_RUNS}
                runs["xfoil"] = run_xfoil(workdir, re, 0.0, TRIPPED_ALPHAS, [],
                                          x.shape[1], UIUC_TRIP)
                for key, (points, _) in runs.items():
                    pts = [points[a] for a in TRIPPED_ALPHAS]
                    out[key] = dict(points=pts)
                    out["summary"][key] = dict(
                        cd_error(measured, pts),
                        converged=sum(p["converged"] for p in pts))
                    s = out["summary"][key]
                    print(f"{out['name']:16s} {key:18s} converged "
                          f"{s['converged']:2d}/{len(pts)}  Cd error mean "
                          f"{s['mean']} rms {s['rms']} ({s['points']} points)",
                          flush=True)
                cases.append(out)
    return cases


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
