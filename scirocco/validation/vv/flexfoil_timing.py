#!/usr/bin/env python3
"""Viscous solve times of Scirocco and flexfoil on the same panel nodes.

    flexfoil_timing.py TIMER DRIVER

TIMER is flexfoil_timing built from flexfoil_timing.cpp and DRIVER is
vv_polar (see README.md). mfoil panels each airfoil; Scirocco takes the
nodes unchanged and flexfoil builds its panels from them. flexfoil runs in
this interpreter, single-threaded except for the batch polar. Prints the
numbers on the documentation's Performance page.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "mfoil"))
import mfoil  # noqa: E402
from flexfoil import _rustfoil as ff  # noqa: E402

REPS = 11
# (naca, panels, Re, alpha)
POINTS = [
    ("0012", 159, 1e6, 2.0), ("0012", 159, 1e6, 8.0),
    ("0012", 199, 1e6, 2.0), ("0012", 199, 1e6, 5.0), ("0012", 199, 1e6, 8.0),
    ("0012", 283, 1e6, 5.0),
    ("2412", 199, 3e6, 4.0), ("4412", 199, 2e5, 4.0),
]
POLAR_ALPHAS = [float(a) for a in range(0, 16)]


def nodes(naca, npanel, workdir):
    m = mfoil.mfoil(naca=naca, npanel=npanel)
    path = workdir / f"nodes_{naca}_{npanel}.txt"
    np.savetxt(path, np.column_stack((m.foil.x.T, m.foil.s, m.foil.t.T)),
               fmt="%.17g")
    # flexfoil takes a flat list from the upper trailing edge, as XFOIL.
    return path, [float(v) for v in m.foil.x.T[::-1].ravel()]


def median_ms(f, reps):
    times, result = [], None
    for _ in range(reps):
        t0 = time.perf_counter()
        result = f()
        times.append(1e3 * (time.perf_counter() - t0))
    return float(np.median(times)), result


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("timer", type=Path)
    parser.add_argument("driver", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        print("| Case | Panels | flexfoil | Scirocco cold | Scirocco warm |")
        for naca, npanel, re, alpha in POINTS:
            path, flat = nodes(naca, npanel, workdir)
            solve = lambda: ff.analyze_faithful(flat, alpha, reynolds=re,
                                                ncrit=9.0, max_iterations=100)
            solve()  # load the library before timing
            t_ff, r = median_ms(solve, REPS)
            out = subprocess.run(
                [str(args.timer), str(path), str(re), str(alpha), str(REPS)],
                capture_output=True, text=True, check=True).stdout.split()
            cold, warm, steps = float(out[0]), float(out[1]), int(out[2])
            conv = r["converged"] and out[3] == "1"
            dcl = abs(float(out[4]) - r["cl"]) / abs(r["cl"])
            print(f"| NACA {naca}, Re {re:.0e}, α {alpha:g}° | {npanel} "
                  f"| {t_ff:.0f} ({r['iterations']}) | {cold:.0f} ({steps}) "
                  f"| {warm:.1f} |{'' if conv else ' NOT CONVERGED'} "
                  f"Cl differs {100 * dcl:.2f}%")

        path, flat = nodes("0012", 199, workdir)
        a = POLAR_ALPHAS
        t, r = median_ms(lambda: [ff.analyze_faithful(flat, x, reynolds=1e6)
                                  for x in a], 5)
        print("\nPolar, NACA 0012, 199 panels, Re 1e6, α 0–15°:")
        print(f"flexfoil one by one   {t:.0f} ms, "
              f"{sum(p['converged'] for p in r)}/{len(a)} converged")
        t, r = median_ms(lambda: ff.analyze_faithful_batch(flat, a, reynolds=1e6), 5)
        print(f"flexfoil batch        {t:.0f} ms on {os.cpu_count()} cores, "
              f"{sum(p['converged'] for p in r)}/{len(a)} converged")
        out = workdir / "polar.json"
        cmd = [str(args.driver), str(path), "1e6", "0", "9",
               ",".join(map(str, a)), "-", str(out)]
        t, _ = median_ms(lambda: subprocess.run(cmd, check=True), 5)
        pts = json.loads(out.read_text())["points"]
        print(f"Scirocco sweep        {t:.0f} ms including process start, "
              f"{sum(p['converged'] for p in pts)}/{len(a)} converged")


if __name__ == "__main__":
    main()
