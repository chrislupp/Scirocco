#!/usr/bin/env python3
"""Produce CSV, figures, and an HTML report from compare_solvers.py outputs."""
import argparse
import csv
import html
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SOLVERS = ("scirocco", "vanilla", "xfoil")
LABELS = dict(scirocco="Scirocco", vanilla="Vanilla mfoil", xfoil="XFOIL 6.99")
COLORS = dict(scirocco="#1266ad", vanilla="#dc7e20", xfoil="#9239a7")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    root = args.directory
    records = sorted(
        json.loads((root / "results.json").read_text()), key=lambda r: r["name"]
    )
    by_name = {r["name"]: r for r in records}
    rows = []
    for r in records:
        for solver in SOLVERS:
            value = r[solver]
            row = {k: r[k] for k in ("name", "naca", "panels", "alpha", "re", "ma")}
            row.update(solver=solver, converged=value["converged"])
            row.update(
                {
                    k: value.get(k)
                    for k in ("cl", "cd", "cm", "xtr_upper", "xtr_lower", "residual")
                }
            )
            rows.append(row)
    write_csv(root / "coefficients.csv", rows)
    errors = []
    for r in records:
        if not r["scirocco"]["converged"]:
            continue
        sc = np.genfromtxt(
            root / r["name"] / "scirocco_surface.csv", delimiter=",", names=True
        )
        for solver in ("vanilla", "xfoil"):
            if not r[solver]["converged"]:
                continue
            other = np.genfromtxt(
                root / r["name"] / (solver + "_surface.csv"), delimiter=",", names=True
            )
            row = dict(name=r["name"], reference=solver)
            for key in ("cl", "cd", "cm", "xtr_upper", "xtr_lower"):
                row["delta_" + key] = r["scirocco"].get(key, 0) - r[solver].get(key, 0)
            for key in ("cp", "cf", "theta", "dstar"):
                if key not in other.dtype.names:
                    row[key + "_rms"] = row[key + "_max"] = None
                    continue
                delta = sc[key] - other[key]
                row[key + "_rms"] = float(np.sqrt(np.mean(delta**2)))
                row[key + "_max"] = float(np.max(abs(delta)))
            errors.append(row)
    write_csv(root / "differences.csv", errors)
    statistics = {
        "case_count": len(records),
        "converged": {s: sum(r[s]["converged"] for r in records) for s in SOLVERS},
    }
    for solver in ("vanilla", "xfoil"):
        subset = [r for r in errors if r["reference"] == solver]
        statistics[solver] = {}
        for q in ("cl", "cd", "cm"):
            worst = max(subset, key=lambda r: abs(r["delta_" + q]))
            statistics[solver][q] = dict(
                max_abs=abs(worst["delta_" + q]), case=worst["name"]
            )
    (root / "statistics.json").write_text(json.dumps(statistics, indent=2))

    plt.rcParams.update(
        {"font.size": 10, "axes.grid": True, "grid.alpha": 0.2, "figure.dpi": 140}
    )
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), layout="constrained")
    for i, naca in enumerate(("0012", "2412")):
        for j, key in enumerate(("cl", "cd", "cm")):
            ax = axes[i, j]
            for ma in (0.0, 0.3):
                selected = sorted(
                    [
                        r
                        for r in records
                        if r["naca"] == naca
                        and r["panels"] == 199
                        and r["re"] == 1e6
                        and r["ma"] == ma
                    ],
                    key=lambda r: r["alpha"],
                )
                for solver in SOLVERS:
                    valid = [r for r in selected if r[solver]["converged"]]
                    ax.plot(
                        [r["alpha"] for r in valid],
                        [r[solver][key] for r in valid],
                        color=COLORS[solver],
                        linestyle="-" if ma == 0 else "--",
                        marker="o" if solver == "vanilla" else ".",
                        markerfacecolor="none",
                        label=f"{LABELS[solver]}, M={ma:g}",
                    )
            ax.set(xlabel="Angle [degrees]", ylabel=key.upper(), title=f"NACA {naca}")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Shared 199-panel geometry · Re = 1 million · Ncrit = 9 · cold starts")
    fig.savefig(root / "polars.png")
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), layout="constrained")
    for j, key in enumerate(("cl", "cd", "cm")):
        for solver in SOLVERS:
            selected = [by_name[f"paper_n{n}"] for n in (99, 199, 299)]
            axes[j].plot(
                [r["panels"] for r in selected],
                [r[solver][key] for r in selected],
                "o-" if solver != "vanilla" else "s",
                color=COLORS[solver],
                markerfacecolor="none" if solver == "vanilla" else COLORS[solver],
                label=LABELS[solver],
            )
        axes[j].set(xlabel="Panel count", ylabel=key.upper())
    axes[0].legend()
    fig.suptitle("Panel refinement: NACA 2412 · alpha 2° · Re = 1 million · Mach 0.4")
    fig.savefig(root / "refinement.png")
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), layout="constrained")
    for solver in SOLVERS:
        data = np.genfromtxt(
            root / "paper_n199" / (solver + "_surface.csv"), delimiter=",", names=True
        )
        le = int(np.argmin(data["x"]))
        for j, key in enumerate(("cp", "cf", "theta")):
            for sl, style in ((slice(None, le + 1), "--"), (slice(le, None), "-")):
                axes[j].plot(
                    data["x"][sl],
                    data[key][sl],
                    style,
                    color=COLORS[solver],
                    linewidth=1.4 if solver != "vanilla" else 0.8,
                    label=LABELS[solver] if style == "-" else None,
                )
            axes[j].set(
                xlabel="x/c",
                ylabel={"cp": "Cp", "cf": "Cf", "theta": "Momentum thickness"}[key],
            )
    axes[0].invert_yaxis()
    axes[0].legend()
    fig.suptitle("NACA 2412 benchmark distributions · solid upper / dashed lower")
    fig.savefig(root / "distributions.png")
    plt.close(fig)

    selected = ["0012_a3_m0.3", "paper_n199", "paper_n299", "2412_low_re", "2412_a10"]
    coefficient_rows = []
    for name in selected:
        r = by_name[name]
        for solver in SOLVERS:
            value = r[solver]
            coefficient_rows.append(
                [
                    name,
                    LABELS[solver],
                    str(value["converged"]),
                    *[
                        f"{value[k]:.9f}" if value["converged"] else "—"
                        for k in ("cl", "cd", "cm")
                    ],
                ]
            )
    text = f"""<!doctype html><html><head><meta charset="utf-8"><title>Three-solver comparison</title>
<style>body{{max-width:1150px;margin:40px auto;font:16px system-ui;line-height:1.5;color:#172333}}table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums;font-size:14px}}td,th{{padding:7px;border-bottom:1px solid #ddd;text-align:right}}td:first-child,th:first-child{{text-align:left}}img{{width:100%}}code{{background:#eee;padding:2px}}h2{{margin-top:35px}}</style></head><body>
<h1>Scirocco, vanilla mfoil, and XFOIL</h1><p>2026-09-22 · {len(records)} cases · unit chord · shared exact panel nodes · quarter-chord moment reference.</p>
<p>Converged: Scirocco {statistics['converged']['scirocco']}/{len(records)}, vanilla mfoil {statistics['converged']['vanilla']}/{len(records)}, XFOIL {statistics['converged']['xfoil']}/{len(records)}.</p>
<h2>Protocol and provenance</h2><p>Vanilla mfoil is the untouched 2023-06-28 Python source archived at <code>4ff9e63:temp/mfoil.py</code>, previously matched to <a href="https://websites.umich.edu/~kfid/codes/mfoil/mfoil.py">Fidkowski's published source</a>. NumPy 1.26.4 and SciPy 1.11.4 avoid its legacy sparse-comparison incompatibility. No source or monkey-patch fixes are applied. Only public run settings change: plotting off, Ncrit 9, and maximum iterations 200.</p>
<p>Scirocco uses the vanilla geometry, arc labels, and tangents directly; XFOIL 6.99 uses LOAD/PCOP, with no repaneling. Binary dump coordinates are checked against the supplied nodes to 2e-15. Both surfaces use natural transition (XFOIL trips at x/c=1). Every main case starts cold. Original stopping criteria remain: mfoil residual tolerance 1e-10, Scirocco 1e-11 with its default PTC, and stock XFOIL convergence. Convergence measures are not identical between codes.</p>
<p>XFOIL values come from its full-precision unformatted polar dump, using the PLXINI/PLXADD layout in the <a href="https://web.mit.edu/drela/Public/web/xfoil/xfoil6.99.tgz">official 6.99 source</a>. Text polars are retained and cross-checked; coefficient comparisons use the full-precision values. The unchanged vanilla source and XFOIL polar-format source are retained in <code>sources/</code>; provenance, raw logs, common coordinates, and per-case results are also retained with this report.</p>
<h2>Convergence findings</h2><p>Of the 30 cold-start cases, Scirocco fails at NACA 2412, Re=1e6, Mach=0.4, alpha=1.995° and 1.9975°, with 199 panels. Both converge when its pseudo-transient continuation (PTC) is disabled, reaching residual norms 1.15e-13 and 4.96e-13 respectively. These diagnostic reruns are preserved in <code>diagnostics/</code> and do not replace the failures in the main sweep. Vanilla mfoil and XFOIL converge at both angles. This identifies a Scirocco convergence robustness gap under its default settings.</p>
<h2>Representative coefficients</h2>{table(['Case','Solver','Converged','Cl','Cd','Cm'],coefficient_rows)}
<p>Unconverged values are excluded from agreement statistics. The NACA 2412 / 10° vanilla solve reaches 200 iterations with residual approximately 20; its unconverged forces are not an alternative physical prediction.</p>
<h2>Maximum coefficient differences over converged pairs</h2>{table(['Comparison','Quantity','Maximum absolute difference','Case'],[[LABELS[s]+' vs Scirocco',q.upper(),f"{statistics[s][q]['max_abs']:.9g}",statistics[s][q]['case']] for s in ('vanilla','xfoil') for q in ('cl','cd','cm')])}
<img src="polars.png" alt="Three solver coefficient comparisons">
<h2>Panel refinement and distributions</h2><p>The 99-, 199-, and 299-panel runs use separately generated common meshes. These three levels show resolution sensitivity; they are not a proof of asymptotic mesh convergence. <a href="https://websites.umich.edu/~kfid/codes/mfoil/mfoil.pdf">Fidkowski's paper</a> documents intentional differences from XFOIL, including stagnation treatment and amplification modeling, and reports Cl=0.4889, Cd=0.00617, Cm=-0.0501 for this mfoil benchmark. The published XFOIL mesh is not our common mesh.</p>
<img src="refinement.png" alt="Panel refinement"><img src="distributions.png" alt="Surface distributions">
<p><strong>Confirmed reporting defect:</strong> Scirocco's stored transition x coordinates remain zero. Its actual transition locations in these comparisons are reconstructed by evaluating the transition residual on a copy of the converged state and interpolating between the adjacent surface nodes. Stored and reconstructed values are both preserved in the per-case JSON. No production solver changes were made.</p>
<p>Data: <a href="coefficients.csv">coefficients.csv</a>, <a href="differences.csv">differences.csv</a>, <a href="results.json">results.json</a>, <a href="provenance.json">provenance.json</a>.</p></body></html>"""
    (root / "report.html").write_text(text)
    print(json.dumps(statistics, indent=2))
    print(root / "report.html")


def table(headers, rows):
    return (
        "<table><thead><tr>"
        + "".join("<th>" + html.escape(str(v)) + "</th>" for v in headers)
        + "</tr></thead><tbody>"
        + "".join(
            "<tr>"
            + "".join("<td>" + html.escape(str(v)) + "</td>" for v in row)
            + "</tr>"
            for row in rows
        )
        + "</tbody></table>"
    )


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
