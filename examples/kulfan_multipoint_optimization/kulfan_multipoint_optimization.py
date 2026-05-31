"""Kulfan/CST multipoint viscous airfoil optimization (OpenMDAO + Scirocco).

Builds on examples/kulfan_viscous_optimization but optimizes a single CST
airfoil across **multiple operating points** simultaneously. The geometry
discipline is shared; each design point gets its own analysis discipline
instance so the (alpha, Re, Ma) state and the converged BL solution don't
collide.

Problem statement
-----------------
Minimize a weighted sum of Cd evaluated at N design points,

    f(cst) = Σ_i w_i · Cd_i(cst, alpha_i)

subject to
    Cl_i(cst, alpha_i) = Cl_target_i      ∀ i
    area(cst)          ≥ 0.8 · area_init
    |A_{i+2} − 2 A_{i+1} + A_i| ≤ SMOOTHNESS_MAX  per surface

Architecture
------------
    [Geometry server] ──x──┬─→ [Analysis @ point 1] ─→ Cl_1, Cd_1, Cm_1
                           │
                           ├─→ [Analysis @ point 2] ─→ Cl_2, Cd_2, Cm_2
                           │
                           ⋮
                           └─→ [AirfoilArea] ─→ area

Each analysis component has its own RemoteExplicitComponent connected to
its own scirocco_server process. SLSQP sees a vector alpha DV (length N)
and N independent Cl equality constraints alongside the shared geometry
DVs.

Usage
-----
    python kulfan_multipoint_optimization.py
"""

import subprocess
import time
from pathlib import Path

import grpc
import numpy as np
import openmdao.api as om
from philote_mdo.openmdao import RemoteExplicitComponent


# ---------------------------------------------------------------------------
# Design points — edit to taste. Weights are renormalized to sum to 1.
# ---------------------------------------------------------------------------

DESIGN_POINTS = [
    # name        target_cl   Re      Ma   weight
    ("low_cl",       0.30,   1.0e6,  0.0,    1.0),
    ("mid_cl",       0.50,   1.0e6,  0.0,    1.5),
    ("high_cl",      0.80,   1.0e6,  0.0,    1.0),
]
N_POINTS = len(DESIGN_POINTS)

# ---------------------------------------------------------------------------
# Geometry configuration — must match KulfanDiscipline option defaults
# ---------------------------------------------------------------------------

NPANEL = 100
N_CST_UPPER = 8
N_CST_LOWER = 8
NACA_TMAX = 0.12
TE_OFFSET = 0.5 * NACA_TMAX * (
    0.2969 - 0.126 - 0.3516 + 0.2843 - 0.1015
) / 0.2 * 2.0
TE_UPPER = +TE_OFFSET / 2.0
TE_LOWER = -TE_OFFSET / 2.0

SMOOTHNESS_MAX = 0.05

# ---------------------------------------------------------------------------
# gRPC addresses + binary paths
# ---------------------------------------------------------------------------

GEOMETRY_ADDR = "localhost:50052"
# Analysis ports start at 50051 and skip 50052 (reserved for geometry).
def _analysis_addr(i):
    base = 50051 if i == 0 else 50052 + i
    return f"localhost:{base}"

ANALYSIS_BINARY = "scirocco/build-philote/src/philote/scirocco_server"
GEOMETRY_BINARY = "scirocco/build-philote/src/philote/scirocco_kulfan_server"


# ---------------------------------------------------------------------------
# CST helpers (mirror the inviscid example)
# ---------------------------------------------------------------------------

def _binom(n, i):
    from math import comb
    return float(comb(n, i))


def _bernstein_matrix(psi, n_cst):
    psi = np.asarray(psi)
    order = n_cst - 1
    B = np.zeros((psi.size, n_cst))
    for i in range(n_cst):
        B[:, i] = (
            _binom(order, i) * np.power(psi, i)
            * np.power(1.0 - psi, order - i)
        )
    return B


def _class_function(psi):
    psi = np.asarray(psi)
    return np.sqrt(psi) * (1.0 - psi)


def _naca_half_thickness(x, tmax=NACA_TMAX):
    return (tmax / 0.2) * (
        0.2969 * np.sqrt(x)
        - 0.1260 * x
        - 0.3516 * x ** 2
        + 0.2843 * x ** 3
        - 0.1015 * x ** 4
    )


def _fit_cst_to_naca(n_cst, te_offset, tmax=NACA_TMAX, n_samples=200):
    psi = np.linspace(0.0, 1.0, n_samples)
    B = _bernstein_matrix(psi, n_cst)
    C = _class_function(psi)[:, None]
    A_mat = C * B
    rhs = _naca_half_thickness(psi, tmax) - psi * te_offset
    cst, *_ = np.linalg.lstsq(A_mat, rhs, rcond=None)
    return cst


def _initial_cps():
    cst_u = _fit_cst_to_naca(N_CST_UPPER, TE_UPPER)
    cst_l = _fit_cst_to_naca(N_CST_LOWER, -TE_LOWER)
    return cst_u, cst_l


# ---------------------------------------------------------------------------
# Airfoil area component (same as the single-point examples)
# ---------------------------------------------------------------------------

class AirfoilArea(om.ExplicitComponent):
    def initialize(self):
        self.options.declare("n_nodes", default=NPANEL + 1, types=int)

    def setup(self):
        nn = self.options["n_nodes"]
        self.add_input("x", shape=(2 * nn,), units="m")
        self.add_output("area", val=0.0, units="m**2")
        self.declare_partials("area", "x")

    def compute(self, inputs, outputs):
        nn = self.options["n_nodes"]
        xc = inputs["x"][:nn]
        zc = inputs["x"][nn:]
        xn = np.roll(xc, -1)
        zn = np.roll(zc, -1)
        outputs["area"] = 0.5 * np.abs(np.sum(xc * zn - xn * zc))

    def compute_partials(self, inputs, partials):
        nn = self.options["n_nodes"]
        xc = inputs["x"][:nn]
        zc = inputs["x"][nn:]
        xn = np.roll(xc, -1)
        zn = np.roll(zc, -1)
        cross_sum = np.sum(xc * zn - xn * zc)
        sign = np.sign(cross_sum) if cross_sum != 0.0 else 1.0
        jac = np.zeros(2 * nn)
        jac[:nn] = 0.5 * sign * (zn - np.roll(zc, 1))
        jac[nn:] = 0.5 * sign * (np.roll(xc, 1) - xn)
        partials["area", "x"] = jac


# ---------------------------------------------------------------------------
# Server management
# ---------------------------------------------------------------------------

def _start_server(binary, address):
    proc = subprocess.Popen(
        [binary, address], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    time.sleep(1.5)
    if proc.poll() is not None:
        _, err = proc.communicate()
        raise RuntimeError(
            f"{binary} ({address}) failed to start: {err.decode()}"
        )
    return proc


def start_servers():
    """Launch one geometry server + N analysis servers (one per design point)."""
    geometry = _start_server(GEOMETRY_BINARY, GEOMETRY_ADDR)
    analyses = []
    for i in range(N_POINTS):
        analyses.append(_start_server(ANALYSIS_BINARY, _analysis_addr(i)))
    time.sleep(1.0)
    return geometry, analyses


# ---------------------------------------------------------------------------
# Problem construction
# ---------------------------------------------------------------------------

def _build_problem(geometry_channel, analysis_channels, hdf5_path=None):
    """Wire one geometry + N analyses + area + aggregated objective."""
    nn = NPANEL + 1

    geometry_kwargs = dict(
        channel=geometry_channel,
        n_cst_upper=N_CST_UPPER,
        n_cst_lower=N_CST_LOWER,
        npanel=NPANEL,
        te_upper=TE_UPPER,
        te_lower=TE_LOWER,
    )

    prob = om.Problem()
    model = prob.model

    model.add_subsystem(
        "geometry",
        RemoteExplicitComponent(**geometry_kwargs),
        promotes_outputs=["x", "smoothness_upper", "smoothness_lower"],
        promotes_inputs=["cst_upper", "cst_lower"],
    )
    model.add_subsystem(
        "area_comp",
        AirfoilArea(n_nodes=nn),
        promotes_inputs=["x"],
    )

    # One analysis component per design point, each connected to its own
    # server. The shared geometry's x feeds them all.
    for i, ((name, _cl, _re, _ma, _w), channel) in enumerate(
        zip(DESIGN_POINTS, analysis_channels)
    ):
        analysis_kwargs = dict(
            channel=channel,
            naca="0012",
            npanel=NPANEL,
        )
        if hdf5_path is not None:
            analysis_kwargs["hdf5_file"] = (
                f"{hdf5_path.with_suffix('')}_{name}.h5"
            )
        model.add_subsystem(
            f"analysis_{i}",
            RemoteExplicitComponent(**analysis_kwargs),
            promotes_inputs=["x"],
        )

    # Aggregated drag objective:  cd_total = Σ w_i · cd_i  (weights normalized).
    # Pass weights as ExecComp input defaults so they're constants, not DVs.
    weights = np.array([p[4] for p in DESIGN_POINTS], dtype=float)
    weights /= weights.sum()
    expr_terms = [f"w_{i} * cd_{i}" for i in range(N_POINTS)]
    expr = "cd_total = " + " + ".join(expr_terms)
    exec_kwargs = {"cd_total": 0.0}
    for i in range(N_POINTS):
        exec_kwargs[f"cd_{i}"] = 0.0
        exec_kwargs[f"w_{i}"] = float(weights[i])
    cd_aggregator = om.ExecComp(expr, **exec_kwargs)
    model.add_subsystem("cd_aggregator", cd_aggregator)
    for i in range(N_POINTS):
        model.connect(f"analysis_{i}.cd", f"cd_aggregator.cd_{i}")

    return prob, weights


def _set_baseline_inputs(prob, cst_u, cst_l, alphas):
    prob.set_val("cst_upper", cst_u)
    prob.set_val("cst_lower", cst_l)
    for i, (_name, _cl, re, ma, _w) in enumerate(DESIGN_POINTS):
        prob.set_val(f"analysis_{i}.alpha", float(alphas[i]))
        prob.set_val(f"analysis_{i}.Re", float(re))
        prob.set_val(f"analysis_{i}.Ma", float(ma))


# ---------------------------------------------------------------------------
# Examples
# ---------------------------------------------------------------------------

def example_single_point_check(geometry_channel, analysis_channels):
    """Sanity check: evaluate the NACA-0012-equivalent shape at every point."""
    prob, weights = _build_problem(geometry_channel, analysis_channels)
    prob.setup()
    cst_u, cst_l = _initial_cps()
    alphas = [5.0] * N_POINTS
    _set_baseline_inputs(prob, cst_u, cst_l, alphas)
    prob.run_model()

    print("=" * 72)
    print("Multipoint baseline (NACA-0012-equivalent CST, alpha = 5° at each)")
    print("=" * 72)
    for i, (name, target_cl, re, _ma, w) in enumerate(DESIGN_POINTS):
        cl = prob.get_val(f"analysis_{i}.cl")[0]
        cd = prob.get_val(f"analysis_{i}.cd")[0]
        cm = prob.get_val(f"analysis_{i}.cm")[0]
        print(f"  {name:>10s}  Cl_target={target_cl:.2f}  Re={re:.0e}  "
              f"w={w:.2f}   Cl={cl:.4f}  Cd={cd:.5f}  Cm={cm:+.4f}")
    print(f"  weights (normalized) = {weights}")
    print(f"  cd_total = {prob.get_val('cd_aggregator.cd_total')[0]:.6f}")
    print(f"  area     = {prob.get_val('area_comp.area')[0]:.6f}")
    print("=" * 72)
    return prob


def _reference_area(geometry_channel, analysis_channels, cst_u, cst_l):
    prob, _ = _build_problem(geometry_channel, analysis_channels)
    prob.setup()
    _set_baseline_inputs(prob, cst_u, cst_l, [0.0] * N_POINTS)
    prob.run_model()
    return (
        float(prob.get_val("area_comp.area")[0]),
        prob.get_val("x").copy(),
    )


def example_optimization(geometry_channel, analysis_channels):
    """Minimize weighted Cd at all design points, Cl equality per point."""
    nn = NPANEL + 1
    hdf5_path = Path(__file__).with_name(
        "kulfan_multipoint_optimization_history.h5"
    )

    cst_u_init, cst_l_init = _initial_cps()
    init_area, x_init = _reference_area(
        geometry_channel, analysis_channels, cst_u_init, cst_l_init,
    )

    prob, weights = _build_problem(
        geometry_channel, analysis_channels, hdf5_path=hdf5_path,
    )

    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options["optimizer"] = "SLSQP"
    prob.driver.options["tol"] = 1e-6
    prob.driver.options["maxiter"] = 200
    prob.driver.options["disp"] = True

    prob.model.add_design_var(
        "cst_upper", lower=0.0, upper=0.40, ref=0.20,
    )
    prob.model.add_design_var(
        "cst_lower", lower=0.0, upper=0.40, ref=0.20,
    )
    # One alpha per design point.
    for i, (_name, target_cl, _re, _ma, _w) in enumerate(DESIGN_POINTS):
        prob.model.add_design_var(
            f"analysis_{i}.alpha", lower=-2.0, upper=10.0, ref=5.0,
        )

    prob.model.add_objective("cd_aggregator.cd_total", ref=0.01)
    for i, (_name, target_cl, _re, _ma, _w) in enumerate(DESIGN_POINTS):
        prob.model.add_constraint(
            f"analysis_{i}.cl", equals=target_cl, ref=target_cl,
        )

    min_area = 0.8 * init_area
    prob.model.add_constraint(
        "area_comp.area", lower=min_area, ref=init_area,
    )
    prob.model.add_constraint(
        "smoothness_upper", lower=-SMOOTHNESS_MAX, upper=SMOOTHNESS_MAX,
        ref=SMOOTHNESS_MAX,
    )
    prob.model.add_constraint(
        "smoothness_lower", lower=-SMOOTHNESS_MAX, upper=SMOOTHNESS_MAX,
        ref=SMOOTHNESS_MAX,
    )

    prob.setup()
    _set_baseline_inputs(prob, cst_u_init, cst_l_init,
                         alphas=[5.0] * N_POINTS)

    print("\n" + "=" * 72)
    print(f"Multipoint optimization: minimize Σ w_i·Cd_i across {N_POINTS} "
          f"points")
    for i, (name, target_cl, re, _ma, _w) in enumerate(DESIGN_POINTS):
        print(f"  point {i} ({name:>10s}): Cl = {target_cl:.2f} at Re = "
              f"{re:.0e}, weight = {weights[i]:.3f}")
    print(f"  Area constraint: area >= {min_area:.6f} "
          f"(80% of NACA-0012 init = {init_area:.6f})")
    print(f"  Smoothness constraint: |A_{{i+2}}-2A_{{i+1}}+A_i| <= "
          f"{SMOOTHNESS_MAX} per surface")
    print("=" * 72 + "\n")

    prob.run_driver()

    cst_u_opt = prob.get_val("cst_upper").copy()
    cst_l_opt = prob.get_val("cst_lower").copy()
    alphas = [float(prob.get_val(f"analysis_{i}.alpha")[0])
              for i in range(N_POINTS)]
    cls = [float(prob.get_val(f"analysis_{i}.cl")[0])
           for i in range(N_POINTS)]
    cds = [float(prob.get_val(f"analysis_{i}.cd")[0])
           for i in range(N_POINTS)]
    cms = [float(prob.get_val(f"analysis_{i}.cm")[0])
           for i in range(N_POINTS)]
    cd_total = float(prob.get_val("cd_aggregator.cd_total")[0])
    area = float(prob.get_val("area_comp.area")[0])
    x_opt = prob.get_val("x").copy()

    print("\n" + "=" * 72)
    print("Optimization Results")
    print("=" * 72)
    print(f"  cst_upper = {np.array2string(cst_u_opt, precision=4)}")
    print(f"  cst_lower = {np.array2string(cst_l_opt, precision=4)}")
    for i, (name, target_cl, re, _ma, _w) in enumerate(DESIGN_POINTS):
        print(f"  point {i} ({name:>10s}, target Cl={target_cl:.2f}): "
              f"alpha={alphas[i]:+.3f}°  Cl={cls[i]:.4f}  Cd={cds[i]:.5f}  "
              f"Cm={cms[i]:+.4f}")
    print(f"  Σ w·Cd  = {cd_total:.6f}")
    print(f"  Area    = {area:.6f} (min = {min_area:.6f})")
    print(f"  Iteration history: {hdf5_path}")
    print("=" * 72)

    plot_results(x_init, x_opt, nn, cst_u_opt, cst_l_opt, cls, cds, alphas)
    return prob


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_results(x_init, x_opt, nn, cst_u, cst_l, cls, cds, alphas):
    import matplotlib.pyplot as plt

    fig, (ax_foil, ax_polar, ax_cps) = plt.subplots(
        3, 1, figsize=(10, 9), gridspec_kw={"height_ratios": [2, 2, 1.5]},
    )

    xi_init, zi_init = x_init[:nn], x_init[nn:]
    xo, zo = x_opt[:nn], x_opt[nn:]
    ax_foil.plot(xi_init, zi_init, "b--", lw=1.2, label="Initial (NACA-0012)")
    ax_foil.plot(xo, zo, "r-", lw=1.8, label="Optimized (multipoint)")
    ax_foil.fill_between(xo, zo, alpha=0.08, color="red")
    ax_foil.set_xlabel("x/c")
    ax_foil.set_ylabel("z/c")
    ax_foil.set_title(
        f"Multipoint CST Optimization across {N_POINTS} design points"
    )
    ax_foil.set_aspect("equal")
    ax_foil.legend(loc="upper right")
    ax_foil.grid(True, alpha=0.3)

    # Polar with the design points highlighted.
    point_names = [p[0] for p in DESIGN_POINTS]
    target_cls = [p[1] for p in DESIGN_POINTS]
    ax_polar.scatter(cds, cls, c="r", s=80, zorder=3, label="Optimized")
    for cl_t, cd, cl, name in zip(target_cls, cds, cls, point_names):
        ax_polar.annotate(name, (cd, cl), xytext=(6, 0),
                          textcoords="offset points", fontsize=9)
        ax_polar.plot([cd, cd], [cl_t, cl], "k:", lw=0.5)
    ax_polar.set_xlabel("Cd")
    ax_polar.set_ylabel("Cl")
    ax_polar.set_title("Optimized drag polar at design points")
    ax_polar.grid(True, alpha=0.3)
    ax_polar.legend(loc="upper left")

    idx_u = np.arange(len(cst_u))
    idx_l = np.arange(len(cst_l))
    ax_cps.plot(idx_u, cst_u, "o-", color="C0", label="cst_upper")
    ax_cps.plot(idx_l, cst_l, "s-", color="C3", label="cst_lower")
    ax_cps.axhline(0.0, color="k", lw=0.5, alpha=0.5)
    ax_cps.set_xlabel("Bernstein coefficient index")
    ax_cps.set_ylabel("CST coefficient value")
    ax_cps.legend(loc="upper right")
    ax_cps.grid(True, alpha=0.3)

    fig.tight_layout()
    out = Path(__file__).with_name("kulfan_multipoint_optimization_result.png")
    fig.savefig(out, dpi=150)
    plt.show()
    print(f"Plot saved to {out}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    geometry_proc, analysis_procs = start_servers()

    try:
        geometry_channel = grpc.insecure_channel(GEOMETRY_ADDR)
        analysis_channels = [
            grpc.insecure_channel(_analysis_addr(i)) for i in range(N_POINTS)
        ]

        example_single_point_check(geometry_channel, analysis_channels)
        example_optimization(geometry_channel, analysis_channels)

        geometry_channel.close()
        for ch in analysis_channels:
            ch.close()
    finally:
        for proc in (geometry_proc, *analysis_procs):
            proc.terminate()
            proc.wait()
