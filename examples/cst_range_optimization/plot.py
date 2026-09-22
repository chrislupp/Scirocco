"""Plot a range benchmark result and optional independent validation polars."""
import argparse
import json
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from geometry import coordinates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument(
        "--baseline",
        type=Path,
        help="Earliest stage summary, for its initial geometry and forces",
    )
    parser.add_argument("--validation", type=Path, nargs="*", default=[])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-unconverged", action="store_true")
    args = parser.parse_args()
    result = json.loads(args.summary.read_text())
    baseline = json.loads(args.baseline.read_text()) if args.baseline else result
    if not result["success"] and not args.allow_unconverged:
        raise RuntimeError("Plot requires a converged optimization")
    candidate_label = "Optimized" if result["success"] else "Unconverged candidate"
    initial = baseline["initial"]["design"]
    final = result["design"]
    gap = result["gap"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), layout="constrained")
    for design, label, color in [
        (initial, "Initial symmetric airfoil", "#777777"),
        (final, f"{candidate_label} airfoil", "#1261a0"),
    ]:
        xy = coordinates(
            design["thickness_coefficients"], design["camber_coefficients"], 1600, gap
        ).reshape(2, -1)
        for ax in axes[0]:
            ax.plot(*xy, label=label, color=color, lw=1.8)
    axes[0, 0].set(title="Continuous geometry", xlabel="x/c", ylabel="y/c")
    axes[0, 1].set(
        title="Leading edge",
        xlabel="x/c",
        ylabel="y/c",
        xlim=(-0.003, 0.07),
        ylim=(-0.05, 0.06),
    )
    for ax in axes[0]:
        ax.set_aspect("equal", adjustable="box")
        ax.legend(fontsize=8)
    for record, label, color, style in [
        (baseline["initial"]["evaluation"], "Initial, design points", "#777777", "o--"),
        (result["evaluation"], f"{candidate_label}, design points", "#1261a0", "o-"),
    ]:
        axes[1, 0].plot(
            record["aero.cl"], record["aero.cd"], style, label=label, color=color, ms=4
        )
        axes[1, 1].plot(
            record["aero.cl"], record["aero.cm"], style, label=label, color=color, ms=4
        )
    colors = iter(
        ["#e68613", "#4c956c", "#b14a8a", "#713f9b", "#627d98", "#a65a2e", "#245c45"]
    )
    for path in args.validation:
        validation = json.loads(path.read_text())
        groups = sorted(
            set((r["panels"], r["Re"], r["ncrit"]) for r in validation["cases"])
        )
        for n, re, ncrit in groups:
            rows = sorted(
                [
                    r
                    for r in validation["cases"]
                    if (r["panels"], r["Re"], r["ncrit"]) == (n, re, ncrit)
                ],
                key=lambda r: r["target"],
            )
            color = next(colors, "#888888")
            label = f'{"Initial" if validation.get("initial_design") else candidate_label}: {n} panels, Re={re/1e6:g}M, Ncrit={ncrit}'
            for ax, key in [(axes[1, 0], "cd"), (axes[1, 1], "cm")]:
                ax.plot(
                    [r["target"] for r in rows],
                    [r[key] if r["success"] else np.nan for r in rows],
                    "-",
                    color=color,
                    lw=1.2,
                    label=label,
                )
    axes[1, 0].set(title="Lift-matched drag", xlabel="$C_l$", ylabel="$C_d$")
    axes[1, 1].set(
        title="Quarter-chord pitching moment", xlabel="$C_l$", ylabel="$C_m$"
    )
    axes[1, 1].axhline(
        -0.08, color="black", lw=0.8, ls=":", label="Design moment limit"
    )
    for ax in axes[1]:
        ax.axvspan(0.3, 0.8, alpha=0.035, color="black")
        ax.legend(fontsize=7, loc="best")
    for ax in axes.flat:
        ax.grid(alpha=0.2)
    fig.suptitle("CST range benchmark — Re = 10⁶, Mach = 0, Ncrit = 9")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for extension in ["png", "pdf"]:
        fig.savefig(args.output.with_suffix("." + extension), dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
