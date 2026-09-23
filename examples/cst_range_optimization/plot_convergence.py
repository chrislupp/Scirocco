"""Plot separate optimality and feasibility histories from a saved IPOPT run."""

import argparse
import csv
import json
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np


def history(summary, directory):
    result = summary["optimizer_result"]
    if "iterations" not in result or "inf_du" not in result["iterations"]:
        raise ValueError("Convergence plots require an IPOPT iteration history")
    logs = list((directory / "cst_range_outputs").glob("*/ipopt*"))
    if len(logs) != 1:
        raise ValueError(f"Expected one IPOPT log below {directory}; found {len(logs)}")
    log = logs[0].read_text()
    rows = []
    for line in log.splitlines():
        match = re.match(
            r"^\s*(\d+r?)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+",
            line,
        )
        if match:
            iteration, _, primal, _ = match.groups()
            if iteration.endswith("r"):
                raise ValueError(
                    "Restoration-phase histories are not supported by this plotter"
                )
            rows.append((int(iteration), float(primal)))
    iterations = np.array([r[0] for r in rows])
    optimality = np.asarray(result["iterations"]["inf_du"])
    if not np.array_equal(iterations, np.arange(len(optimality))):
        raise ValueError("IPOPT log and saved dual-infeasibility history do not align")
    # CasADi's callback inf_pr includes the internal slack formulation. The
    # printed inf_pr instead measures violations of the original NLP constraints.
    feasibility = np.array([r[1] for r in rows])
    footer = re.search(r"Constraint violation\.+:\s+([\d.eE+-]+)\s+([\d.eE+-]+)", log)
    if footer:
        feasibility[-1] = float(footer.group(2))
    for residual in (optimality, feasibility):
        if not np.all(np.isfinite(residual)) or np.any(residual < 0):
            raise ValueError("Residual histories must be finite and nonnegative")
    return iterations, optimality, feasibility


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument(
        "--output", type=Path, help="Output directory; defaults to the run directory"
    )
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text())
    iterations, optimality, feasibility = history(summary, args.summary.parent)
    out = args.output or args.summary.parent
    out.mkdir(parents=True, exist_ok=True)
    with (out / "convergence.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["iteration", "dual_infeasibility", "primal_infeasibility"])
        writer.writerows(zip(iterations, optimality, feasibility))
    for name, values, title, ylabel, color in [
        ("optimality", optimality, "Optimality", "Dual infeasibility", "#0072B2"),
        ("feasibility", feasibility, "Feasibility", "Primal infeasibility", "#D55E00"),
    ]:
        fig, ax = plt.subplots(figsize=(8, 4.6), layout="constrained")
        # Zero is preserved in the CSV; mark it at a labeled floor on the log plot.
        positive = values[values > 0]
        floor = min(1e-16, float(np.min(positive)) / 10) if len(positive) else 1e-16
        ax.semilogy(
            iterations, np.maximum(values, floor), "-o", color=color, lw=1.7, ms=3.2
        )
        if np.any(values == 0):
            ax.plot(
                iterations[values == 0],
                np.full(np.sum(values == 0), floor),
                "v",
                color=color,
                label=f"Zero residual shown at {floor:.0e}",
            )
        tolerance = summary["settings"]["optimizer_tol"]
        ax.axhline(
            tolerance,
            color="#666666",
            ls="--",
            lw=1.2,
            label=f"Requested NLP tolerance: {tolerance:g}",
        )
        ax.set(
            title=f"{title} over iterations",
            xlabel="Iteration",
            ylabel=f"{ylabel} (log scale)",
            xlim=(0, max(1, iterations[-1]) + 1),
        )
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid(which="major", alpha=0.22)
        ax.legend(
            loc="lower left" if name == "optimality" else "upper right", fontsize=9
        )
        ax.text(
            0.97 if name == "optimality" else 0.03,
            0.95 if name == "optimality" else 0.15,
            f"Final: {values[-1]:.3e}",
            transform=ax.transAxes,
            ha="right" if name == "optimality" else "left",
            va="top",
            color=color,
        )
        for extension in ["png", "pdf"]:
            fig.savefig(out / f"{name}.{extension}", dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    main()
