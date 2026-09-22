"""Range-based Scirocco benchmark. All output and server processes belong to --output."""
import os

os.environ["OPENMDAO_REPORTS"] = "false"
os.environ["MPLBACKEND"] = "Agg"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import signal
import shutil
import time
import traceback
import numpy as np
import openmdao.api as om
from modopt import IPOPT, SLSQP, OpenMDAOProblem
from optimality import check_kkt
from geometry import CamberThicknessCST, area_weights, bernstein, coordinates, elevate
from parallel_aero import AeroPool, ParallelAero

ROOT = Path(__file__).resolve().parents[2]
TARGETS = np.array([0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85])
WEIGHTS = np.r_[
    0.0, 0.08, 0.20, 0.24, 0.21333333333333333, 0.18666666666666668, 0.08, 0.0
]
WRT = ["thickness_coefficients", "camber_coefficients", "aero.alpha"]
BASE_TARGETS = TARGETS.copy()
BASE_WEIGHTS = WEIGHTS.copy()
SOURCE_FILES = ["run.py", "geometry.py", "parallel_aero.py", "optimality.py"]


def serialize(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    return str(x)


def reference():
    # Same continuous NACA 0012 CST fit used by the existing multipoint example.
    import importlib.util

    path = (
        ROOT
        / "examples/kulfan_multipoint_optimization/kulfan_multipoint_optimization.py"
    )
    spec = importlib.util.spec_from_file_location("original_multipoint", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    upper, lower = module._initial_cps()
    assert np.allclose(upper, lower)
    return upper, float(module.TE_UPPER - module.TE_LOWER)


def match_lift(pool, x, alpha):
    alpha = np.array(alpha, copy=True, dtype=float)
    previous = None
    for iteration in range(15):
        values = pool.evaluate(x, alpha)
        error = values["cl"] - TARGETS
        if np.max(abs(error)) < 1e-9:
            return alpha, {k: v.copy() for k, v in values.items()}
        slope = np.full(len(alpha), 0.1)
        if previous is not None:
            change = alpha - previous[0]
            np.divide(error - previous[1], change, out=slope, where=abs(change) > 1e-12)
        if np.any(~np.isfinite(slope)) or np.any(abs(slope) < 1e-5):
            raise RuntimeError("Invalid lift-matching slope")
        previous = (alpha.copy(), error.copy())
        alpha -= np.clip(error / slope, -1.0, 1.0)
    raise RuntimeError("Could not match baseline lift to 1e-9")


def capture(prob):
    return {
        "design": {
            "thickness_coefficients": prob.get_val(WRT[0]).copy(),
            "camber_coefficients": prob.get_val(WRT[1]).copy(),
            "alpha": prob.get_val(WRT[2]).copy(),
        },
        "evaluation": {
            k: prob.get_val(k).copy()
            for k in [
                "x",
                "geometry.area",
                "geometry.thickness",
                "geometry.wedge",
                "aero.cl",
                "aero.cd",
                "aero.cm",
                "objective.cd_total",
            ]
        },
    }


def derivative_check(prob, pool):
    outputs = ["objective.cd_total", "aero.cl", "aero.cm"]
    base = [prob.get_val(k).copy() for k in WRT]
    rng = np.random.default_rng(7301)
    direction = [rng.normal(size=len(a)) * s for a, s in zip(base, [0.02, 0.02, 0.1])]
    totals = prob.compute_totals(of=outputs, wrt=WRT, return_format="array")
    analytic = totals @ np.concatenate(direction)
    records = []
    try:
        for step in [1e-3, 1e-4, 1e-5]:
            samples = []
            for sign in [1, -1]:
                for k, a, d in zip(WRT, base, direction):
                    prob.set_val(k, a + sign * step * d)
                prob.run_model()
                samples.append(
                    np.concatenate([prob.get_val(k).reshape(-1) for k in outputs])
                )
            finite = (samples[0] - samples[1]) / (2 * step)
            records.append(
                {
                    "step": step,
                    "max_absolute_error": float(np.max(abs(finite - analytic))),
                    "finite_difference": finite.tolist(),
                }
            )
    finally:
        for k, a in zip(WRT, base):
            prob.set_val(k, a)
        prob.run_model()
    passed = sum(r["max_absolute_error"] < 2e-6 for r in records) >= 2
    return {
        "passed": passed,
        "analytic": analytic.tolist(),
        "records": records,
        "description": "One deterministic mixed shape/angle direction, all eight Cl and Cm outputs and weighted Cd; central differences at three steps.",
    }


def main():
    global TARGETS, WEIGHTS

    def stop(signum, frame):
        raise KeyboardInterrupt(f"Received signal {signum}")

    signal.signal(signal.SIGTERM, stop)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--panels", type=int, default=100)
    parser.add_argument("--thickness-degree", type=int, default=4)
    parser.add_argument("--camber-degree", type=int, default=3)
    parser.add_argument("--seed", type=Path)
    parser.add_argument(
        "--allow-unconverged-seed",
        action="store_true",
        help="Explicitly restart from a saved candidate; its status is recorded",
    )
    parser.add_argument("--max-iter", type=int, default=160)
    parser.add_argument("--optimizer", choices=["ipopt", "slsqp"], default="ipopt")
    parser.add_argument("--hessian-history", type=int, default=20)
    parser.add_argument(
        "--start", choices=["symmetric", "cambered", "thicker"], default="symmetric"
    )
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument(
        "--extra-lift-point",
        type=float,
        action="append",
        default=[],
        help="Add a zero-objective-weight constraint point within the guard interval",
    )
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--skip-check", action="store_true")
    parser.add_argument(
        "--binary",
        type=Path,
        default=Path("/tmp/scirocco-gradient-philote/src/philote/scirocco_server"),
    )
    args = parser.parse_args()
    if any(not 0.25 <= point <= 0.85 for point in args.extra_lift_point):
        parser.error("Extra lift points must lie within [0.25, 0.85]")
    TARGETS = np.unique(np.r_[BASE_TARGETS, args.extra_lift_point])
    WEIGHTS = np.zeros(len(TARGETS))
    for point, weight in zip(BASE_TARGETS, BASE_WEIGHTS):
        WEIGHTS[np.argmin(abs(TARGETS - point))] = weight
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    args.binary = args.binary.resolve()
    if args.seed:
        args.seed = args.seed.resolve()
    if (out / "summary.json").exists():
        raise RuntimeError(
            "Choose a fresh output directory; results are never overwritten"
        )
    os.chdir(out)
    (out / "source").mkdir(exist_ok=True)
    for filename in SOURCE_FILES:
        shutil.copy2(Path(__file__).parent / filename, out / "source" / filename)
    started = time.monotonic()
    pool = None
    prob = None
    optimizer = None
    summary = {
        "success": False,
        "commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "settings": {
            **vars(args),
            "targets": TARGETS,
            "weights": WEIGHTS,
            "Re": 1e6,
            "Mach": 0,
            "ncrit": 9,
            "optimizer_tol": 1e-6,
            "radius_bounds": [0.005, 0.02],
            "moment_lower": -0.08,
            "area_fraction": 0.8,
            "thickness_fraction": 0.8,
            "wedge_bounds": [4, 20],
        },
        "binary_sha256": hashlib.sha256(args.binary.read_bytes()).hexdigest(),
        "source_sha256": {
            f: hashlib.sha256((Path(__file__).parent / f).read_bytes()).hexdigest()
            for f in SOURCE_FILES
        },
    }
    try:
        ref, gap = reference()
        summary["reference_thickness_coefficients"] = ref
        summary["gap"] = gap
        p, q = args.thickness_degree, args.camber_degree
        fit_x = np.linspace(0, 1, 1001)
        c = np.sqrt(fit_x) * (1 - fit_x)
        t = np.linalg.lstsq(
            c[:, None] * bernstein(fit_x, p),
            c * (bernstein(fit_x, len(ref) - 1) @ ref),
            rcond=None,
        )[0]
        m = np.zeros(q + 1)
        if args.start == "cambered":
            t *= 1 + 0.03 * np.cos(np.linspace(0, np.pi, len(t)))
            m[:] = 0.02
        elif args.start == "thicker":
            t *= 1.03
        alpha = TARGETS / 0.11
        reference_alpha_guess = alpha.copy()
        if args.seed:
            seed = json.loads(args.seed.resolve().read_text())
            summary["seed_optimizer_success"] = bool(seed["success"])
            if not seed["success"] and not args.allow_unconverged_seed:
                raise RuntimeError(
                    "A continuation seed must be a converged optimization"
                )
            t = elevate(seed["design"]["thickness_coefficients"], p)
            m = elevate(seed["design"]["camber_coefficients"], q)
            alpha = np.interp(
                TARGETS, seed["settings"]["targets"], seed["design"]["alpha"]
            )
            reference_alpha_guess = np.interp(
                TARGETS,
                seed["settings"]["targets"],
                seed["reference_aerodynamics"]["alpha"],
            )
        pool = AeroPool(args.binary, args.panels, TARGETS, out, args.workers)
        refalpha, refvalues = match_lift(
            pool, coordinates(ref, np.zeros(1), args.panels, gap), reference_alpha_guess
        )
        summary["reference_aerodynamics"] = {"alpha": refalpha, **refvalues}
        alpha, initial_values = match_lift(
            pool, coordinates(t, m, args.panels, gap), alpha
        )
        area_ref = float(area_weights(len(ref) - 1) @ ref + gap / 2)
        locations = np.r_[np.linspace(0.15, 0.65, 21), 0.9]
        thickness_ref = (
            2
            * np.sqrt(locations)
            * (1 - locations)
            * (bernstein(locations, len(ref) - 1) @ ref)
            + gap * locations
        )
        prob = om.Problem()
        prob.model.add_subsystem(
            "geometry",
            CamberThicknessCST(
                panels=args.panels,
                thickness_degree=p,
                camber_degree=q,
                gap=gap,
                thickness_stations=locations,
            ),
            promotes_inputs=WRT[:2],
            promotes_outputs=["x"],
        )
        prob.model.add_subsystem("aero", ParallelAero(pool=pool), promotes_inputs=["x"])
        prob.model.add_subsystem(
            "objective",
            om.ExecComp(
                "cd_total=dot(weights,cd)", weights=WEIGHTS, cd=np.zeros(len(TARGETS))
            ),
        )
        prob.model.connect("aero.cd", "objective.cd")
        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options["disp"] = False
        prob.model.add_design_var(
            WRT[0],
            lower=np.r_[0.1, np.zeros(p)],
            upper=np.r_[0.2, np.full(p, 0.4)],
            ref=0.2,
        )
        prob.model.add_design_var(WRT[1], lower=-0.4, upper=0.4, ref=0.1)
        prob.model.add_design_var(WRT[2], lower=-2, upper=10, ref=5)
        prob.model.add_objective("objective.cd_total", ref=0.01)
        prob.model.add_constraint("aero.cl", equals=TARGETS, ref=TARGETS)
        prob.model.add_constraint("aero.cm", lower=-0.08, ref=0.1)
        prob.model.add_constraint(
            "aero.cd",
            indices=[0, len(TARGETS) - 1],
            upper=1.05 * refvalues["cd"][[0, -1]],
            ref=0.01,
        )
        # Keep these in the total constraint Jacobian: the installed modopt
        # adapter omits OpenMDAO constraints marked linear=True from its totals.
        prob.model.add_constraint("geometry.area", lower=0.8 * area_ref, ref=area_ref)
        prob.model.add_constraint(
            "geometry.thickness", lower=0.8 * thickness_ref, ref=thickness_ref
        )
        prob.model.add_constraint("geometry.wedge", lower=4, upper=20, ref=10)
        prob.setup()
        for key, value in zip(WRT, [t, m, alpha]):
            prob.set_val(key, value)
        prob.run_model()
        summary["initial"] = capture(prob)
        summary["reference_area"] = area_ref
        summary["reference_thickness"] = thickness_ref
        summary["thickness_stations"] = locations
        print(
            "INITIAL",
            json.dumps(
                {
                    "design": summary["initial"]["design"],
                    "cd": initial_values["cd"],
                    "weighted_cd": summary["initial"]["evaluation"][
                        "objective.cd_total"
                    ],
                },
                default=serialize,
            ),
            flush=True,
        )
        if not args.skip_check:
            summary["derivative_check"] = derivative_check(prob, pool)
            (out / "derivative_check.json").write_text(
                json.dumps(summary["derivative_check"], indent=2) + "\n"
            )
            print(
                "DERIVATIVE CHECK", json.dumps(summary["derivative_check"]), flush=True
            )
            if not summary["derivative_check"]["passed"]:
                raise RuntimeError("Aerodynamic directional derivative check failed")
        if args.check_only:
            summary["checks_complete"] = True
            return
        model = OpenMDAOProblem(problem_name="cst_range", om_problem=prob)
        if args.optimizer == "ipopt":
            optimizer = IPOPT(
                model,
                solver_options={
                    "max_iter": args.max_iter,
                    "tol": 1e-6,
                    "print_level": 5,
                    "hessian_approximation": "limited-memory",
                    "limited_memory_max_history": args.hessian_history,
                },
            )
        else:
            optimizer = SLSQP(
                model,
                solver_options={"maxiter": args.max_iter, "ftol": 1e-12, "disp": True},
                readable_outputs=["x"],
            )
        result = optimizer.solve()
        summary["optimizer_result"] = result
        model.check_if_warm_and_run_model(result["x"], force_rerun=True)
        summary.update(capture(prob))
        summary["kkt_check"] = check_kkt(model, result["x"])
        summary["success"] = bool(result["success"]) and summary["kkt_check"]["passed"]
        print(
            "RESULT",
            json.dumps(
                {
                    k: v
                    for k, v in summary.items()
                    if k in ["success", "design", "evaluation"]
                },
                default=serialize,
            ),
            flush=True,
        )
    except BaseException:
        summary["error"] = traceback.format_exc()
        print(summary["error"], flush=True)
        if optimizer is not None and hasattr(optimizer, "results"):
            summary["optimizer_result"] = optimizer.results
        if prob is not None:
            try:
                summary["last_evaluation"] = capture(prob)
            except Exception:
                pass
    finally:
        summary["total_seconds"] = time.monotonic() - started
        if pool is not None:
            summary["analysis_batches"] = pool.calls
            summary["derivative_batches"] = pool.derivative_calls
            pool.close()
        if prob is not None:
            prob.cleanup()
        (out / "summary.json").write_text(
            json.dumps(summary, indent=2, default=serialize) + "\n"
        )
    if not summary["success"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
