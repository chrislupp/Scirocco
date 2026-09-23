#!/usr/bin/env python3
"""Run unmodified published Python mfoil, Scirocco, and XFOIL on shared nodes.

Run with NumPy 1.26.4 / SciPy 1.11.4, which support vanilla mfoil's legacy
sparse comparison. No monkey patches or source changes are applied to mfoil.
The archived file must be provided explicitly; its hash is recorded.

Example:
  c++ -O3 -std=c++17 -I scirocco/include -isystem /opt/homebrew/include/eigen3 \
    scirocco/validation/compare_solver_driver.cpp /tmp/scirocco-gradient-philote/libscirocco.a \
    -o /tmp/scirocco-compare
  OPENBLAS_NUM_THREADS=1 /tmp/scirocco-vanilla-mfoil-env/bin/python \
    scirocco/validation/compare_solvers.py --vanilla /tmp/scirocco-three-solver/vanilla_mfoil.py \
    --driver /tmp/scirocco-compare --output pa_output/three_solver_2026-09-22
"""
import argparse
import contextlib
import hashlib
import importlib.util
import json
import subprocess
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import scipy
from scipy.io import FortranFile


def load_vanilla(path):
    spec = importlib.util.spec_from_file_location("vanilla_mfoil", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def cases():
    result = []

    def add(name, naca, alpha, re=1e6, ma=0.0, panels=199, gradients=False):
        result.append(
            dict(
                name=name,
                naca=naca,
                alpha=alpha,
                re=re,
                ma=ma,
                panels=panels,
                gradients=gradients,
            )
        )

    for naca in ("0012", "2412"):
        for ma in (0.0, 0.3):
            for alpha in (0.0, 3.0, 6.0):
                add(f"{naca}_a{alpha:g}_m{ma:g}", naca, alpha, ma=ma)
            add(f"{naca}_inviscid_m{ma:g}", naca, 3, re=0, ma=ma)
        add(f"{naca}_a10", naca, 10)
    for panels in (99, 199, 299):
        add(
            f"paper_n{panels}",
            "2412",
            2,
            ma=0.4,
            panels=panels,
        )
    add("2412_low_re", "2412", 4, re=2e5)
    for alpha in (1.98, 1.99, 1.995, 1.9975, 2.0025, 2.005, 2.01, 2.02):
        add(f"paper_alpha_{alpha:g}", "2412", alpha, ma=0.4)
    return result


def run_case(case, options):
    started = time.monotonic()
    directory = Path(options["output"]) / case["name"]
    directory.mkdir(parents=True, exist_ok=True)
    module = load_vanilla(options["vanilla"])
    model = module.mfoil(naca=case["naca"], npanel=case["panels"])
    model.param.verb = 1
    model.param.doplot = False
    model.param.niglob = 200  # public iteration limit, same as other solvers
    model.param.ncrit = 9
    # Keep vanilla's own panels, spline arc labels, and tangents exactly.
    geometry = np.column_stack((model.foil.x.T, model.foil.s, model.foil.t.T))
    np.savetxt(directory / "geometry.txt", geometry, fmt="%.17g")
    # XFOIL conventional ordering is upper TE -> LE -> lower TE.
    np.savetxt(
        directory / "airfoil.dat",
        model.foil.x.T[::-1],
        fmt="%.17g",
        header=f"common_NACA{case['naca']}",
        comments="",
    )
    record = dict(case)
    try:
        with (directory / "vanilla.log").open("w") as log, contextlib.redirect_stdout(
            log
        ), contextlib.redirect_stderr(log):
            model.setoper(alpha=case["alpha"], Ma=case["ma"])
            if case["re"]:
                model.setoper(Re=case["re"])
            model.solve()
            v = dict(
                converged=bool(model.glob.conv),
                cl=float(model.post.cl),
                cd=float(model.post.cd),
                cm=float(model.post.cm),
                cdf=float(model.post.cdf),
                cdpi=float(model.post.cdpi),
                rtol=model.param.rtol,
                residual=float(np.linalg.norm(model.glob.R)) if case["re"] else 0.0,
            )
            if case["re"]:
                v.update(
                    xtr_upper=float(model.vsol.Xt[1, 1]),
                    xtr_lower=float(model.vsol.Xt[0, 1]),
                )
                module.get_distributions(model)
            n = model.foil.N
            surface = np.column_stack(
                (
                    model.foil.x.T,
                    model.post.cp[:n],
                    model.post.cf[:n] if case["re"] else np.zeros(n),
                    model.post.th[:n] if case["re"] else np.zeros(n),
                    model.post.ds[:n] if case["re"] else np.zeros(n),
                )
            )
            np.savetxt(
                directory / "vanilla_surface.csv",
                surface,
                delimiter=",",
                header="x,y,cp,cf,theta,dstar",
                comments="",
            )
        record["vanilla"] = v
    except Exception:
        record["vanilla"] = dict(converged=False, error=traceback.format_exc())
    command = [
        options["driver"],
        str(directory / "geometry.txt"),
        str(case["alpha"]),
        str(case["re"]),
        str(case["ma"]),
        str(directory),
        str(int(case["gradients"])),
    ]
    with (directory / "scirocco.log").open("w") as log:
        try:
            completed = subprocess.run(command, stdout=log, stderr=log, timeout=180)
            assert completed.returncode == 0, f"driver exit {completed.returncode}"
            record["scirocco"] = json.loads((directory / "scirocco.json").read_text())
        except Exception:
            record["scirocco"] = dict(converged=False, error=traceback.format_exc())
    commands = [
        "PLOP",
        "G F",
        "",
        "LOAD airfoil.dat",
        "PCOP",
        "PSAV xfoil_nodes.dat",
        "OPER",
        f"MACH {case['ma']}",
        "ITER 200",
        "VPAR",
        "N 9",
        "XTR 1 1",
        "",
    ]
    if case["re"]:
        commands.append(f"VISC {case['re']}")
    commands.extend(
        [
            "PACC",
            "polar.txt",
            "polar_dump.bin",
            f"ALFA {case['alpha']}",
            "CPWR cp.txt",
            "DUMP bl.txt",
            "",
            "QUIT",
            "",
        ]
    )
    script = "\n".join(commands)
    (directory / "xfoil.in").write_text(script)
    # The XFOIL file-writing commands prompt before overwriting existing files.
    for name in ("polar.txt", "polar_dump.bin", "cp.txt", "bl.txt", "xfoil_nodes.dat"):
        (directory / name).unlink(missing_ok=True)
    try:
        xf = subprocess.run(
            [options["xfoil"]],
            input=script,
            text=True,
            capture_output=True,
            cwd=directory,
            timeout=90,
        )
        (directory / "xfoil.log").write_text(xf.stdout + xf.stderr)
        lines = (directory / "polar.txt").read_text().splitlines()
        start = next(
            i for i, line in enumerate(lines) if line.strip().startswith("------")
        )
        rows = [
            list(map(float, line.split()))
            for line in lines[start + 1 :]
            if line.strip()
        ]
        if not rows:
            raise RuntimeError("No converged polar row")
        row = rows[-1]
        assert abs(row[0] - case["alpha"]) < 5.1e-4
        assert xf.returncode == 0, f"XFOIL exit {xf.returncode}"
        # Stock XFOIL PLXINI / PLXADD (src/xpol.f) writes double-precision
        # forces and BL distributions to an unformatted Fortran polar dump.
        # Read that instead of differentiating rounded text polar columns.
        with FortranFile(directory / "polar_dump.bin", "r") as dump:
            header = dump.read_record(np.uint8)
            assert bytes(header[48:56]) == b"XFOIL   "
            dump.read_reals(np.float64)  # Mach, Re/1e6, Ncrit(upper,lower)
            dump.read_ints(np.int32)  # Mach/Re variation types
            info = dump.read_ints(np.int32)
            assert info[3] == model.foil.N
            xfnodes = dump.read_reals(np.float64).reshape(-1, 2)
            exact = dump.read_reals(np.float64)
            assert len(exact) == 7
            assert abs(exact[0] - case["alpha"]) < 1e-12
            for j, tolerance in ((1, 5.1e-5), (2, 5.1e-6), (4, 5.1e-5)):
                assert abs(exact[j] - row[j]) < tolerance
            bl_info = dump.read_ints(np.int32)
            upper = dump.read_reals(np.float64).reshape(-1, 6)
            lower = dump.read_reals(np.float64).reshape(-1, 6)
        assert xfnodes.shape == model.foil.x.T.shape
        error = np.max(abs(xfnodes - model.foil.x.T[::-1]))
        assert error < 2e-15, f"XFOIL changed geometry: {error}"
        cp = np.loadtxt(directory / "cp.txt", comments="#")
        assert len(cp) == model.foil.N
        surface = np.column_stack((model.foil.x.T, cp[::-1, -1]))
        header = "x,y,cp"
        if case["re"]:
            # Each side starts with a synthetic stagnation state (index 0),
            # followed by the original surface nodes and then wake points.
            bl = np.vstack((upper[1 : bl_info[2]][::-1], lower[1 : bl_info[3]]))
            assert len(bl) == model.foil.N
            np.testing.assert_allclose(bl[:, 0], xfnodes[:, 0], atol=2e-15, rtol=0)
            np.testing.assert_allclose(bl[:, 1], cp[:, -1], atol=5.1e-6, rtol=0)
            surface = np.column_stack(
                (model.foil.x.T, bl[::-1, 1], bl[::-1, 4], bl[::-1, 2], bl[::-1, 3])
            )
            header = "x,y,cp,cf,theta,dstar"
        np.savetxt(
            directory / "xfoil_surface.csv",
            surface,
            delimiter=",",
            header=header,
            comments="",
        )
        record["xfoil"] = dict(
            converged=True,
            cl=exact[1],
            cd=exact[2],
            cm=exact[4],
            xtr_upper=exact[5],
            xtr_lower=exact[6],
            saved_coordinate_error=float(error),
        )
    except Exception:
        record["xfoil"] = dict(converged=False, error=traceback.format_exc())
    record["seconds"] = time.monotonic() - started
    (directory / "results.json").write_text(json.dumps(record, indent=2))
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vanilla", required=True, type=Path)
    parser.add_argument("--driver", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--xfoil", default="/opt/homebrew/bin/xfoil")
    parser.add_argument("--workers", default=3, type=int)
    parser.add_argument("--only", help="comma-separated case names")
    args = parser.parse_args()
    options = vars(args).copy()
    for key in ("vanilla", "driver", "output"):
        options[key] = str(options[key].resolve())
    args.output.mkdir(parents=True, exist_ok=True)
    provenance = dict(
        options,
        vanilla_sha256=hashlib.sha256(args.vanilla.read_bytes()).hexdigest(),
        numpy=np.__version__,
        scipy=scipy.__version__,
        scirocco_commit=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        vanilla_version="2023-06-28",
        vanilla_git_source="4ff9e63:temp/mfoil.py",
        vanilla_url="https://websites.umich.edu/~kfid/codes/mfoil/mfoil.py",
        xfoil_version="6.99",
        moment_reference=[0.25, 0],
        chord=1,
        ncrit=9,
        xfoil_forced_transition=[1, 1],
        cold_start=True,
        max_iterations=200,
        scirocco_ptc=True,
        xfoil_repanel=False,
    )
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2))
    selected = [
        c for c in cases() if not args.only or c["name"] in args.only.split(",")
    ]
    records = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_case, c, options): c for c in selected}
        for future in as_completed(futures):
            record = future.result()
            records.append(record)
            print(
                record["name"],
                {key: record[key] for key in ("scirocco", "vanilla", "xfoil")},
                flush=True,
            )
            (args.output / "results.json").write_text(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
