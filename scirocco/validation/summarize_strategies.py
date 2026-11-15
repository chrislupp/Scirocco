#!/usr/bin/env python3
"""Summarize the staged solver experiments and isolated final timing runs."""
import argparse
import html
import json
from pathlib import Path


def table(headers, rows):
    def cells(tag, row):
        return (
            "<tr>"
            + "".join(f"<{tag}>{html.escape(str(x))}</{tag}>" for x in row)
            + "</tr>"
        )

    return (
        "<table>"
        + cells("th", headers)
        + "".join(cells("td", r) for r in rows)
        + "</table>"
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    args = p.parse_args()
    root = args.directory
    names = ["baseline", "scaled", "trust", "topology", "continuation"]
    labels = [
        "Original Newton/PTC",
        "1: equation and variable scaling",
        "2: dogleg, continuous update acceptance",
        "3: dogleg, complete state acceptance",
        "4: scaled Newton/PTC with continuation fallback",
    ]
    values = [json.loads((root / n / "summary.json").read_text()) for n in names]
    rows = [
        [
            label,
            f"{s['original']['converged']}/30",
            f"{s['extended']['converged']}/32",
            f"{s['converged']}/62",
            s["iterations"],
        ]
        for label, s in zip(labels, values)
    ]
    timings = []
    for name, label in [
        ("timing_baseline", "Original"),
        ("timing_final", "Final default"),
    ]:
        if (root / name / "summary.json").exists():
            s = json.loads((root / name / "summary.json").read_text())
            timings.append(
                [label, f"{s['converged']}/62", f"{s['seconds']:.2f}", s["iterations"]]
            )
    records = json.loads(
        (
            root
            / (
                "timing_final"
                if (root / "timing_final/results.json").exists()
                else "continuation"
            )
            / "results.json"
        ).read_text()
    )
    recovered = [
        [
            r["name"],
            r["continuation_attempts"],
            r["continuation_steps"],
            r["continuation_iterations"],
            f"{r['cl']:.10g}",
            f"{r['cd']:.10g}",
            f"{r['cm']:.10g}",
        ]
        for r in records
        if r.get("continuation_attempts")
    ]
    links = ", ".join(f'<a href="{name}/results.csv">{name} CSV</a>' for name in names)
    text = f"""<!doctype html><meta charset="utf-8"><title>Solver strategy benchmarks</title>
<style>body{{max-width:1150px;margin:35px auto;font:16px system-ui;line-height:1.5;color:#172333}}table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}}td,th{{padding:7px;text-align:right;border-bottom:1px solid #ddd}}td:first-child,th:first-child{{text-align:left}}</style>
<h1>Solver strategy benchmarks</h1>
<p>Scaling and adaptive continuation are enabled by default. Trust regions remain opt-in because their full-state implementation converged in fewer cases than the scaled Newton/PTC method. The final strategy converged in all 62 tested cases; this is not an unconditional convergence guarantee.</p>
<h2>Protocol</h2><p>All runs start cold on the identical vanilla-mfoil panel coordinates used in the three-solver comparison. The original 30 cases are supplemented by 32 fixed cases: NACA 0012 and 2412, 199 panels, alpha -4, 8, 12, 15 degrees, Reynolds numbers 200,000 and 1,000,000, Mach 0 and 0.4. These difficult cases test numerical convergence, not physical validity outside the model's range. No XFOIL/mfoil reference solutions were computed for the additional 32 cases.</p>
<p>The physical equations and 1e-11 full coupled residual tolerance are preserved. Continuation keeps the target wake geometry and approaches the exact requested angle or lift target. Its additional work is bounded by 24 attempts and 800 coupled iterations; successful direct solves bypass it. Trial failures retain the last accepted checkpoint. Continuation is currently unavailable when explicit wake rebuilding is enabled.</p>
<h2>Staged results</h2>{table(['Strategy','Original cases','Additional cases','All cases','Residual evaluations in outer loops'], rows)}
<p>The outer-loop count includes the convergence check and every continuation attempt, but excludes internal trust-region trial residual evaluations. It therefore understates trust-region work. Stage 2 checks the continuous update before topology changes; it is retained only for reproducing the staged experiment. Stage 3 includes those changes in acceptance. Both trust-region variants bypass PTC even if its configuration flag is set.</p>
<h2>Isolated timing comparison</h2>{table(['Strategy','Converged','Elapsed seconds','Outer-loop evaluations'], timings)}
<p>Timing runs are sequential, one process at a time, after builds and tests completed. Times include process startup, solve, surface output and transition reconstruction. Each solver is measured once, so small timing differences are not significant. Earlier exploratory timings overlapped build/test activity and are not used for the cost comparison.</p>
<h2>Cases recovered by continuation</h2>{table(['Case','Attempts','Accepted steps','Extra iterations','Cl','Cd','Cm'], recovered)}
<p>Across cases where both scaling-only and continuation-stage runs converged, maximum coefficient differences were Cl={values[-1]['max_coefficient_differences']['cl']:.3g}, Cd={values[-1]['max_coefficient_differences']['cd']:.3g}, Cm={values[-1]['max_coefficient_differences']['cm']:.3g}.</p>
<p>Raw per-case logs, geometries by reference, JSON and surface distributions are in each stage directory. Reproduction tools are in <code>scirocco/validation/</code>. Data: {links}.</p>"""
    (root / "report.html").write_text(text)
    print(root / "report.html")


if __name__ == "__main__":
    main()
