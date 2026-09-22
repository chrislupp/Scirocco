#!/usr/bin/env python3
"""Benchmark one driver on the fixed 30-case comparison plus 32 harder cases.

Use --reference with an earlier output directory to report paired differences.
All cases start cold, use the same nodes, and retain the 1e-11 stopping tolerance.
Run stages sequentially (one worker) for comparable elapsed times.
"""
import argparse
import csv
import hashlib
import json
import re
import subprocess
import time
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--driver", type=Path, required=True)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--reference", type=Path)
    p.add_argument("--strategy", type=int)
    p.add_argument("--ptc", type=int, default=1)
    p.add_argument("--only")
    args = p.parse_args()
    source = args.source.resolve()
    cases = json.loads((source / "results.json").read_text())
    cases = [
        {
            **{k: r[k] for k in ("name", "naca", "panels", "alpha", "re", "ma")},
            "geometry": str(source / r["name"] / "geometry.txt"),
            "suite": "original",
        }
        for r in cases
    ]
    for naca in ("0012", "2412"):
        geometry = str(source / f"{naca}_a0_m0" / "geometry.txt")
        for re_value in (2e5, 1e6):
            for ma in (0.0, 0.4):
                for alpha in (-4.0, 8.0, 12.0, 15.0):
                    cases.append(
                        dict(
                            name=f"extended_{naca}_a{alpha:g}_r{re_value:g}_m{ma:g}",
                            naca=naca,
                            panels=199,
                            alpha=alpha,
                            re=re_value,
                            ma=ma,
                            geometry=geometry,
                            suite="extended",
                        )
                    )
    if args.only:
        cases = [r for r in cases if r["name"] in args.only.split(",")]
    args.output.mkdir(parents=True, exist_ok=True)
    previous = {}
    if args.reference:
        previous = {
            r["name"]: r
            for r in json.loads((args.reference / "results.json").read_text())
        }
    results = []
    for r in cases:
        directory = (args.output / r["name"]).resolve()
        directory.mkdir(exist_ok=True)
        cmd = [
            str(args.driver.resolve()),
            r["geometry"],
            str(r["alpha"]),
            str(r["re"]),
            str(r["ma"]),
            str(directory),
            "0",
            str(args.ptc),
        ]
        if args.strategy is not None:
            cmd.append(str(args.strategy))
        start = time.monotonic()
        with (directory / "scirocco.log").open("w") as log:
            try:
                proc = subprocess.run(cmd, stdout=log, stderr=log, timeout=120)
                value = (
                    json.loads((directory / "scirocco.json").read_text())
                    if proc.returncode == 0
                    else dict(converged=False, exit_code=proc.returncode)
                )
            except subprocess.TimeoutExpired:
                value = dict(converged=False, timeout=True)
        row = {**r, **value, "seconds": time.monotonic() - start}
        log = (directory / "scirocco.log").read_text()
        row["iterations"] = len(re.findall(r"Iteration \d+:", log))
        if r["name"] in previous:
            old = previous[r["name"]]
            row["previous_converged"] = old["converged"]
            if old["converged"] and row["converged"]:
                for k in ("cl", "cd", "cm"):
                    row["delta_" + k] = row[k] - old[k]
        results.append(row)
        (args.output / "results.json").write_text(json.dumps(results, indent=2))
        print(
            row["name"],
            row["converged"],
            row["iterations"],
            round(row["seconds"], 3),
            flush=True,
        )
    summary = dict(
        driver_sha256=hashlib.sha256(args.driver.read_bytes()).hexdigest(),
        strategy=args.strategy,
        ptc=args.ptc,
        cases=len(results),
        converged=sum(r["converged"] for r in results),
        seconds=sum(r["seconds"] for r in results),
        iterations=sum(r["iterations"] for r in results),
        failures=[r["name"] for r in results if not r["converged"]],
    )
    for suite in ("original", "extended"):
        subset = [r for r in results if r["suite"] == suite]
        summary[suite] = dict(
            cases=len(subset), converged=sum(r["converged"] for r in subset)
        )
    if previous:
        summary["recovered"] = [
            r["name"]
            for r in results
            if r["converged"] and not previous[r["name"]]["converged"]
        ]
        summary["regressed"] = [
            r["name"]
            for r in results
            if not r["converged"] and previous[r["name"]]["converged"]
        ]
        summary["max_coefficient_differences"] = {
            k: max(
                [abs(r["delta_" + k]) for r in results if "delta_" + k in r] or [0.0]
            )
            for k in ("cl", "cd", "cm")
        }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2))
    fields = sorted(set().union(*(r.keys() for r in results)))
    with (args.output / "results.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
