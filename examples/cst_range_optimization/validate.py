"""Cold-start lift-matched validation of a saved range benchmark result."""
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from scipy.optimize import minimize_scalar
from geometry import area_weights, bernstein, coordinates


def geometry_checks(summary):
    t = np.array(summary["design"]["thickness_coefficients"])
    m = np.array(summary["design"]["camber_coefficients"])
    ref = np.array(summary["reference_thickness_coefficients"])
    gap = summary["gap"]

    def thickness(x, coeff):
        return (
            2
            * np.sqrt(x)
            * (1 - x)
            * (bernstein(np.atleast_1d(x), len(coeff) - 1) @ coeff)
            + gap * x
        )

    def ratio(x):
        return float((thickness(x, t) / thickness(x, ref))[0])

    # Refine each station interval independently, also testing the endpoints.
    grid = np.linspace(0.15, 0.65, 2001)
    fractions = thickness(grid, t) / thickness(grid, ref)
    candidates = [(float(fractions[i]), float(grid[i])) for i in [0, len(grid) - 1]]
    for a, b in zip(np.linspace(0.15, 0.65, 21)[:-1], np.linspace(0.15, 0.65, 21)[1:]):
        opt = minimize_scalar(
            ratio, bounds=(a, b), method="bounded", options={"xatol": 1e-13}
        )
        candidates.append((float(opt.fun), float(opt.x)))
    candidates.append((float(np.min(fractions)), float(grid[np.argmin(fractions)])))
    minimum, where = min(candidates)
    up, lo = -m[-1] - t[-1] + gap / 2, -m[-1] + t[-1] - gap / 2
    return {
        "radius": float(t[0] ** 2 / 2),
        "area_fraction": float(
            (area_weights(len(t) - 1) @ t + gap / 2) / summary["reference_area"]
        ),
        "minimum_thickness_fraction": minimum,
        "minimum_thickness_x": where,
        "thickness_fraction_at_90_percent": ratio(0.9),
        "wedge_degrees": float(np.degrees(np.arctan(lo) - np.arctan(up))),
        "minimum_thickness_on_dense_grid": float(
            np.min(thickness(np.linspace(0.000001, 1, 10001), t))
        ),
    }


def check_point(driver, geometry, target, guess, re, ncrit, directory, timeout):
    directory.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    history = []
    alpha = float(guess)
    previous = None
    try:
        for iteration in range(12):
            result = directory / f"solve_{iteration}.json"
            with (directory / f"solve_{iteration}.log").open("w") as log:
                subprocess.run(
                    [
                        str(driver),
                        str(geometry),
                        str(alpha),
                        str(re),
                        "0",
                        str(ncrit),
                        str(result),
                    ],
                    stdout=log,
                    stderr=log,
                    check=True,
                    timeout=timeout,
                )
            value = json.loads(result.read_text())
            history.append(value)
            if not value["converged"]:
                raise RuntimeError("Flow did not converge")
            error = value["cl"] - target
            if abs(error) < 1e-8:
                return {
                    "success": True,
                    "target": target,
                    "Re": re,
                    "ncrit": ncrit,
                    **value,
                    "seconds": time.monotonic() - started,
                    "history": history,
                }
            slope = (
                0.1
                if previous is None or abs(alpha - previous[0]) < 1e-12
                else (error - previous[1]) / (alpha - previous[0])
            )
            if not np.isfinite(slope) or abs(slope) < 1e-5:
                raise RuntimeError("Invalid lift-matching slope")
            previous = (alpha, error)
            alpha -= float(np.clip(error / slope, -1.0, 1.0))
        raise RuntimeError("Lift matching exceeded 12 cold solves")
    except Exception as exc:
        return {
            "success": False,
            "target": target,
            "Re": re,
            "ncrit": ncrit,
            "error": str(exc),
            "seconds": time.monotonic() - started,
            "history": history,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--driver", type=Path, default=Path("/tmp/scirocco-cst-validate")
    )
    parser.add_argument("--mode", choices=["mesh", "polar", "stress"], default="mesh")
    parser.add_argument("--panels", type=int, nargs="+", default=[200, 400, 800])
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--allow-unconverged",
        action="store_true",
        help="Check candidate function values without certifying optimizer convergence",
    )
    parser.add_argument(
        "--initial",
        action="store_true",
        help="Validate the initial design instead of the optimum",
    )
    parser.add_argument("--timeout", type=float, default=180)
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text())
    if not args.initial and not summary["success"] and not args.allow_unconverged:
        raise RuntimeError("Validation requires a converged optimization")
    if args.initial:
        summary["design"] = summary["initial"]["design"]
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out / "summary.json").exists():
        raise RuntimeError("Choose a fresh validation output directory")
    driver = args.driver.resolve()
    started = time.monotonic()
    data = {
        "source_summary": str(args.summary.resolve()),
        "driver_sha256": hashlib.sha256(driver.read_bytes()).hexdigest(),
        "mode": args.mode,
        "initial_design": args.initial,
        "source_optimizer_success": summary["success"],
        "geometry": geometry_checks(summary),
        "cases": [],
    }
    targets = (
        np.array(summary["settings"]["targets"])
        if args.mode == "mesh"
        else np.linspace(0.25, 0.85, 25)
    )
    conditions = (
        [(1e6, 9)] if args.mode != "stress" else [(0.7e6, 9), (1.3e6, 9), (1e6, 5)]
    )
    if args.mode == "stress":
        targets = np.array(summary["settings"]["targets"])
    guesses = np.interp(
        targets, summary["settings"]["targets"], summary["design"]["alpha"]
    )
    tasks = {}
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        for panels in args.panels:
            xy = (
                coordinates(
                    summary["design"]["thickness_coefficients"],
                    summary["design"]["camber_coefficients"],
                    panels,
                    summary["gap"],
                )
                .reshape(2, -1)
                .T
            )
            geometry = out / f"geometry_{panels}.dat"
            np.savetxt(geometry, xy, fmt="%.17g")
            for re, ncrit in conditions:
                for i, (target, guess) in enumerate(zip(targets, guesses)):
                    directory = out / f"n{panels}_re{re:g}_nc{ncrit}_cl{target:.3f}"
                    f = executor.submit(
                        check_point,
                        driver,
                        geometry,
                        float(target),
                        guess,
                        re,
                        ncrit,
                        directory,
                        args.timeout,
                    )
                    tasks[f] = panels
        for future in as_completed(tasks):
            record = {"panels": tasks[future], **future.result()}
            data["cases"].append(record)
            print(
                json.dumps({k: v for k, v in record.items() if k != "history"}),
                flush=True,
            )
            data["total_seconds"] = time.monotonic() - started
            data["success"] = all(r["success"] for r in data["cases"]) and len(
                data["cases"]
            ) == len(tasks)
            (out / "summary.json").write_text(json.dumps(data, indent=2) + "\n")
    if not data["success"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
