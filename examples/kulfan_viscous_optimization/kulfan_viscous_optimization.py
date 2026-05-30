"""Kulfan/CST-parameterized viscous airfoil optimization (OpenMDAO + Scirocco).

Minimize Cd subject to a Cl-equality and an area-floor constraint. Same
problem as ``examples/naca_optimization`` (the original 3-DV NACA 4-digit
example), but the geometry is parameterized by Kulfan's CST: a
Bernstein-polynomial shape function times a fixed-exponent class function
(`N1=0.5, N2=1.0`) that builds in a round LE and a closed TE.

Why CST here:
  * The Bernstein basis is non-negative on [0, 1], so positive coefficients
    give a non-negative shape interpolant *everywhere*, not just at sample
    points. Bounding `cst_upper[i] >= 0`, `cst_lower[i] >= 0` is sufficient
    to keep the airfoil physical at every intermediate optimizer iterate.
  * The class function `C(psi) = sqrt(psi) * (1 - psi)` provides the round
    LE automatically (vertical tangent at psi = 0).

Architecture::

    [KulfanDiscipline] --x--> [SiroccoDiscipline] --cl,cd,cm--> [SLSQP]
       (Philote gRPC)  |        (Philote gRPC)
                       +--> [AirfoilArea] --area-->

Two Philote-Cpp gRPC servers run side by side: the geometry server on
50052, the analysis server on 50051. OpenMDAO chains them.

Usage:
    python kulfan_viscous_optimization.py
"""

import subprocess
import time
from pathlib import Path

import grpc
import numpy as np
import openmdao.api as om
from philote_mdo.openmdao import RemoteExplicitComponent


# ---------------------------------------------------------------------------
# Configuration — keep in sync with KulfanDiscipline option defaults
# ---------------------------------------------------------------------------

NPANEL = 100              # panels; nnodes = NPANEL + 1
N_CST_UPPER = 8
N_CST_LOWER = 8
NACA_TMAX = 0.12          # initial NACA-0012-like thickness
TE_OFFSET = 0.5 * NACA_TMAX * (
    0.2969 - 0.126 - 0.3516 + 0.2843 - 0.1015
) / 0.2 * 2.0  # full TE gap from NACA-4-digit formula at x=1
TE_UPPER = +TE_OFFSET / 2.0
TE_LOWER = -TE_OFFSET / 2.0

ANALYSIS_ADDR = "localhost:50051"
GEOMETRY_ADDR = "localhost:50052"

ANALYSIS_BINARY = "scirocco/build-philote/src/philote/scirocco_server"
GEOMETRY_BINARY = "scirocco/build-philote/src/philote/scirocco_kulfan_server"


# ---------------------------------------------------------------------------
# CST helpers (mirror the C++ implementation so we can build initial guesses)
# ---------------------------------------------------------------------------

def _binom(n, i):
    """Binomial coefficient n choose i (as float)."""
    from math import comb
    return float(comb(n, i))


def _bernstein_matrix(psi, n_cst):
    """Bernstein basis matrix [len(psi), n_cst].

    `B[j, i] = K^(n_cst-1)_i * psi_j^i * (1 - psi_j)^(n_cst-1-i)`.
    """
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
    """C(psi) = psi^0.5 * (1 - psi)^1.0 — round LE, sharp TE."""
    psi = np.asarray(psi)
    return np.sqrt(psi) * (1.0 - psi)


def _naca_half_thickness(x, tmax=NACA_TMAX):
    """NACA 4-digit finite-TE half-thickness law."""
    return (tmax / 0.2) * (
        0.2969 * np.sqrt(x)
        - 0.1260 * x
        - 0.3516 * x ** 2
        + 0.2843 * x ** 3
        - 0.1015 * x ** 4
    )


def _fit_cst_to_naca(n_cst, te_offset, tmax=NACA_TMAX, n_samples=200):
    """Least-squares fit CST coefficients to a symmetric NACA half-thickness.

    Solves `(C(psi) * B(psi)) @ A ≈ y_target(psi) - psi * te_offset`.
    """
    psi = np.linspace(0.0, 1.0, n_samples)
    B = _bernstein_matrix(psi, n_cst)
    C = _class_function(psi)[:, None]
    A_mat = C * B
    rhs = _naca_half_thickness(psi, tmax) - psi * te_offset
    cst, *_ = np.linalg.lstsq(A_mat, rhs, rcond=None)
    return cst


def _initial_cps():
    """NACA-0012-equivalent initial guess for upper / lower CST coeffs."""
    cst_u = _fit_cst_to_naca(N_CST_UPPER, TE_UPPER)
    cst_l = _fit_cst_to_naca(N_CST_LOWER, -TE_LOWER)
    return cst_u, cst_l


# ---------------------------------------------------------------------------
# Airfoil area component (shoelace formula, analytical gradients)
# ---------------------------------------------------------------------------

class AirfoilArea(om.ExplicitComponent):
    """Cross-sectional area from flattened ``[x; z]`` panel coordinates."""

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

        dA_dxc = 0.5 * sign * (zn - np.roll(zc, 1))
        dA_dzc = 0.5 * sign * (np.roll(xc, 1) - xn)

        jac = np.zeros(2 * nn)
        jac[:nn] = dA_dxc
        jac[nn:] = dA_dzc
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
    """Launch both analysis and Kulfan-geometry Philote servers."""
    analysis = _start_server(ANALYSIS_BINARY, ANALYSIS_ADDR)
    geometry = _start_server(GEOMETRY_BINARY, GEOMETRY_ADDR)
    time.sleep(1.0)
    return analysis, geometry


# ---------------------------------------------------------------------------
# Problem construction
# ---------------------------------------------------------------------------

def _build_problem(analysis_channel, geometry_channel, hdf5_file=None):
    """Wire KulfanDiscipline → SiroccoDiscipline + area through OpenMDAO."""
    nn = NPANEL + 1

    geometry_kwargs = dict(
        channel=geometry_channel,
        n_cst_upper=N_CST_UPPER,
        n_cst_lower=N_CST_LOWER,
        npanel=NPANEL,
        te_upper=TE_UPPER,
        te_lower=TE_LOWER,
    )
    analysis_kwargs = dict(
        channel=analysis_channel,
        naca="0012",
        npanel=NPANEL,
    )
    if hdf5_file is not None:
        analysis_kwargs["hdf5_file"] = str(hdf5_file)

    prob = om.Problem()
    model = prob.model

    model.add_subsystem(
        "geometry",
        RemoteExplicitComponent(**geometry_kwargs),
        promotes_outputs=["x"],
        promotes_inputs=["cst_upper", "cst_lower"],
    )
    model.add_subsystem(
        "area_comp",
        AirfoilArea(n_nodes=nn),
        promotes_inputs=["x"],
    )
    model.add_subsystem(
        "analysis",
        RemoteExplicitComponent(**analysis_kwargs),
        promotes_inputs=["x"],
    )
    return prob


REYNOLDS = 1e6


def _set_baseline_inputs(prob, cst_u, cst_l, alpha=5.0):
    prob.set_val("cst_upper", cst_u)
    prob.set_val("cst_lower", cst_l)
    prob.set_val("analysis.alpha", alpha)
    prob.set_val("analysis.Re", REYNOLDS)
    prob.set_val("analysis.Ma", 0.0)


# ---------------------------------------------------------------------------
# Examples
# ---------------------------------------------------------------------------

def example_single_point(analysis_channel, geometry_channel):
    """Evaluate the CST geometry tuned to mimic NACA-0012 at alpha = 5°."""
    prob = _build_problem(analysis_channel, geometry_channel)
    prob.setup()

    cst_u, cst_l = _initial_cps()
    _set_baseline_inputs(prob, cst_u, cst_l, alpha=5.0)
    prob.run_model()

    cl = prob.get_val("analysis.cl")[0]
    cd = prob.get_val("analysis.cd")[0]
    cm = prob.get_val("analysis.cm")[0]
    area = prob.get_val("area_comp.area")[0]

    print("=" * 60)
    print("Single-Point Analysis: NACA-0012-equivalent CST shape (viscous)")
    print("=" * 60)
    print(f"  alpha = 5.0 deg, Re = {REYNOLDS:.0e}, Ma = 0.0")
    print(f"  cst_upper = {np.array2string(cst_u, precision=4)}")
    print(f"  cst_lower = {np.array2string(cst_l, precision=4)}")
    print(f"  Cl   = {cl:.6f}")
    print(f"  Cd   = {cd:.6f}")
    print(f"  Cm   = {cm:.6f}")
    print(f"  Area = {area:.6f}")
    print("=" * 60)
    return prob


def _reference_area(analysis_channel, geometry_channel, cst_u, cst_l):
    """Compute area for a given CST configuration."""
    prob = _build_problem(analysis_channel, geometry_channel)
    prob.setup()
    _set_baseline_inputs(prob, cst_u, cst_l, alpha=0.0)
    prob.run_model()
    return (
        float(prob.get_val("area_comp.area")[0]),
        prob.get_val("x").copy(),
    )


def example_optimization(analysis_channel, geometry_channel, target_cl=0.5):
    """Minimize Cd subject to Cl = target_cl and area >= 80% of NACA-0012.

    Viscous (Re = 1e6). Mirrors the original 3-DV NACA optimization but with
    a CST-parameterized geometry (16 CST DVs + alpha).
    """
    nn = NPANEL + 1
    hdf5_path = Path(__file__).with_name(
        "kulfan_viscous_optimization_history.h5"
    )

    cst_u_init, cst_l_init = _initial_cps()
    init_area, x_init = _reference_area(
        analysis_channel, geometry_channel, cst_u_init, cst_l_init,
    )

    prob = _build_problem(analysis_channel, geometry_channel,
                          hdf5_file=hdf5_path)

    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options["optimizer"] = "SLSQP"
    prob.driver.options["tol"] = 1e-6
    prob.driver.options["maxiter"] = 200
    prob.driver.options["disp"] = True

    # Bernstein non-negativity guarantees a physical airfoil whenever all
    # coefficients are >= 0 (with the upper/lower sign convention). The
    # viscous BL solver now uses pseudo-transient continuation, which
    # regularizes the Newton step and keeps the iteration stable even when
    # SLSQP visits geometries with oscillatory surface curvature — so we
    # can use the wider bounds that the parameterization actually permits.
    prob.model.add_design_var(
        "cst_upper", lower=0.0, upper=0.40, ref=0.20,
    )
    prob.model.add_design_var(
        "cst_lower", lower=0.0, upper=0.40, ref=0.20,
    )
    prob.model.add_design_var(
        "analysis.alpha", lower=-2.0, upper=10.0, ref=5.0,
    )

    prob.model.add_objective("analysis.cd", ref=0.01)
    prob.model.add_constraint("analysis.cl", equals=target_cl, ref=target_cl)
    min_area = 0.8 * init_area
    prob.model.add_constraint(
        "area_comp.area", lower=min_area, ref=init_area,
    )

    prob.setup()
    _set_baseline_inputs(prob, cst_u_init, cst_l_init, alpha=5.0)

    print("\n" + "=" * 60)
    print(f"Viscous optimization: minimize Cd, target Cl = {target_cl}, "
          f"Re = {REYNOLDS:.0e}")
    print(f"  DVs: {N_CST_UPPER} upper CST + {N_CST_LOWER} lower CST + alpha")
    print(f"  Area constraint: area >= {min_area:.6f} "
          f"(80% of NACA-0012 init = {init_area:.6f})")
    print("=" * 60 + "\n")

    prob.run_driver()

    cst_u_opt = prob.get_val("cst_upper").copy()
    cst_l_opt = prob.get_val("cst_lower").copy()
    alpha = prob.get_val("analysis.alpha")[0]
    cl = prob.get_val("analysis.cl")[0]
    cd = prob.get_val("analysis.cd")[0]
    cm = prob.get_val("analysis.cm")[0]
    area = prob.get_val("area_comp.area")[0]
    x_opt = prob.get_val("x").copy()

    print("\n" + "=" * 60)
    print("Optimization Results")
    print("=" * 60)
    print(f"  cst_upper = {np.array2string(cst_u_opt, precision=4)}")
    print(f"  cst_lower = {np.array2string(cst_l_opt, precision=4)}")
    print(f"  alpha = {alpha:.4f} deg")
    print(f"  Cl    = {cl:.6f}")
    print(f"  Cd    = {cd:.6f}")
    print(f"  Cm    = {cm:.6f}")
    print(f"  Area  = {area:.6f} (min = {min_area:.6f})")
    print(f"  Iteration history: {hdf5_path}")
    print("=" * 60)

    plot_results(x_init, x_opt, nn, target_cl, cd, alpha,
                 cst_u_opt, cst_l_opt)
    return prob


def plot_results(x_init, x_opt, nn, target_cl, opt_cd, opt_alpha,
                 cst_u, cst_l):
    """Plot initial vs optimized airfoils plus the CST coefficient values."""
    import matplotlib.pyplot as plt

    xi_init, zi_init = x_init[:nn], x_init[nn:]
    xo, zo = x_opt[:nn], x_opt[nn:]

    fig, (ax_foil, ax_cps) = plt.subplots(
        2, 1, figsize=(10, 6), gridspec_kw={"height_ratios": [3, 2]},
    )

    ax_foil.plot(xi_init, zi_init, "b--", lw=1.2, label="Initial (NACA-0012)")
    ax_foil.plot(xo, zo, "r-", lw=1.8, label="Optimized (CST)")
    ax_foil.fill_between(xo, zo, alpha=0.08, color="red")
    ax_foil.set_xlabel("x/c")
    ax_foil.set_ylabel("z/c")
    ax_foil.set_title(
        rf"Viscous CST Optimization — min $C_d$ s.t. $C_l = {target_cl}$ "
        rf"(Re = {REYNOLDS:.0e})"
        "\n"
        rf"Optimized: $C_d = {opt_cd:.5f}$,  "
        rf"$\alpha = {opt_alpha:.2f}°$"
    )
    ax_foil.set_aspect("equal")
    ax_foil.legend(loc="upper right")
    ax_foil.grid(True, alpha=0.3)

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
    out = Path(__file__).with_name("kulfan_viscous_optimization_result.png")
    fig.savefig(out, dpi=150)
    plt.show()
    print(f"Plot saved to {out}")


if __name__ == "__main__":
    analysis_proc, geometry_proc = start_servers()

    try:
        analysis_channel = grpc.insecure_channel(ANALYSIS_ADDR)
        geometry_channel = grpc.insecure_channel(GEOMETRY_ADDR)

        example_single_point(analysis_channel, geometry_channel)
        example_optimization(analysis_channel, geometry_channel, target_cl=0.5)

        analysis_channel.close()
        geometry_channel.close()
    finally:
        for proc in (analysis_proc, geometry_proc):
            proc.terminate()
            proc.wait()
