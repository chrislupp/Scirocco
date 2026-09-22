"""Check every shape and angle derivative at a saved design."""
import os

os.environ["OPENMDAO_REPORTS"] = "false"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import argparse
import json
from pathlib import Path
import time
import numpy as np
import openmdao.api as om
from geometry import CamberThicknessCST
from parallel_aero import AeroPool


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--binary",
        type=Path,
        default=Path("/tmp/scirocco-gradient-philote/src/philote/scirocco_server"),
    )
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text())
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out / "summary.json").exists():
        raise RuntimeError("Choose a fresh output directory")
    d = summary["design"]
    t = np.array(d["thickness_coefficients"])
    m = np.array(d["camber_coefficients"])
    alpha = np.array(d["alpha"])
    prob = om.Problem()
    shape = CamberThicknessCST(
        panels=summary["settings"]["panels"],
        thickness_degree=len(t) - 1,
        camber_degree=len(m) - 1,
        gap=summary["gap"],
    )
    prob.model.add_subsystem("shape", shape, promotes=["*"])
    prob.setup()
    prob.set_val("thickness_coefficients", t)
    prob.set_val("camber_coefficients", m)
    prob.run_model()
    x = prob.get_val("x").copy()
    geometry_jac = np.column_stack([shape.jt, shape.jm])
    pool = AeroPool(
        args.binary.resolve(),
        summary["settings"]["panels"],
        summary["settings"]["targets"],
        out,
    )
    result = {
        "source_summary": str(args.summary.resolve()),
        "source_optimizer_success": summary["success"],
        "records": [],
    }
    started = time.monotonic()
    try:
        jac = pool.derivatives(x, alpha)
        analytic = {
            k: np.column_stack([jac[k, "x"] @ geometry_jac, jac[k, "alpha"]])
            for k in ["cl", "cd", "cm"]
        }
        result["analytic"] = {k: v.tolist() for k, v in analytic.items()}
        result["columns"] = (
            [f"T{i}" for i in range(len(t))]
            + [f"M{i}" for i in range(len(m))]
            + ["own_alpha"]
        )
        # Each force row depends only on its own angle, so a simultaneous
        # angle perturbation checks all eight nonzero angle derivatives.
        for factor in [1.0, 0.1]:
            columns = {k: [] for k in analytic}
            for col in range(geometry_jac.shape[1] + 1):
                angle = col == geometry_jac.shape[1]
                h = factor * (1e-5 if angle else 2e-6)
                samples = []
                for sign in [1, -1]:
                    xp = x if angle else x + sign * h * geometry_jac[:, col]
                    ap = alpha + sign * h if angle else alpha
                    samples.append(pool.evaluate(xp, ap))
                for key in analytic:
                    columns[key].append((samples[0][key] - samples[1][key]) / (2 * h))
            finite = {k: np.column_stack(v) for k, v in columns.items()}
            error = {k: float(np.max(abs(finite[k] - analytic[k]))) for k in analytic}
            record = {
                "step_factor": factor,
                "max_absolute_errors": error,
                "finite_difference": {k: v.tolist() for k, v in finite.items()},
            }
            result["records"].append(record)
            print(
                json.dumps(
                    {k: v for k, v in record.items() if k != "finite_difference"}
                ),
                flush=True,
            )
        # Finite differences are a diagnostic, not a machine-precision oracle.
        result["passed"] = any(
            all(
                r["max_absolute_errors"][k] < limit
                for k, limit in [("cl", 3e-6), ("cd", 3e-8), ("cm", 1e-6)]
            )
            for r in result["records"]
        )
    finally:
        pool.close()
        prob.cleanup()
        result["seconds"] = time.monotonic() - started
        (out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    if not result.get("passed"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
