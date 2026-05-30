"""Spline-parameterized airfoil optimization using OpenMDAO and Scirocco.

Same optimization problem as ``examples/naca_optimization`` — minimize drag
subject to lift and area constraints — but the geometry is parameterized by
a free-form cubic-spline camber line plus a perpendicular thickness
distribution (Scirocco ``GeometryDiscipline``) instead of three NACA digits.

Architecture::

    [GeometryDiscipline] --x--> [SiroccoDiscipline] --cl,cd,cm--> [SLSQP]
        (Philote gRPC)   |        (Philote gRPC)
        n_camber_cp + n_thickness_cp + alpha DVs
                         +--> [AirfoilArea] --area-->

Both disciplines run as separate Philote-Cpp gRPC servers (ports 50051 and
50052). The geometry server provides analytical partials of node
coordinates w.r.t. the spline control points; the analysis server provides
adjoint partials of (cl, cd, cm) w.r.t. node coordinates. OpenMDAO chains
them together.

Usage:
    python spline_optimization.py
"""

import subprocess
import time
from pathlib import Path

import grpc
import numpy as np
import openmdao.api as om
from philote_mdo.openmdao import RemoteExplicitComponent


# ---------------------------------------------------------------------------
# Configuration: must match GeometryDiscipline option defaults / arguments
# ---------------------------------------------------------------------------

NPANEL = 100              # panels per side-and-back; nnodes = NPANEL + 1
N_CAMBER_CP = 6
N_THICKNESS_CP = 8
CLOSED_TE = True

ANALYSIS_ADDR = "localhost:50051"
GEOMETRY_ADDR = "localhost:50052"

ANALYSIS_BINARY = "scirocco/build-philote/src/philote/scirocco_server"
GEOMETRY_BINARY = "scirocco/build-philote/src/philote/scirocco_geometry_server"


# ---------------------------------------------------------------------------
# Helpers — mirror the C++ knot-placement and the NACA thickness law so we
# can build a sensible initial guess and area reference.
# ---------------------------------------------------------------------------

def _xi_thickness_knots(n_thickness_cp=N_THICKNESS_CP, closed_te=CLOSED_TE):
    """Cosine-clustered xi positions of the *full* thickness-spline knot set.

    Mirrors ``GeometryDiscipline::BuildBasisMatrices``: total knots
    = n_thickness_cp + 1 + (closed_te ? 1 : 0), placed cosine-clustered on
    [0, 1]. Free knot k corresponds to index ``1 + k`` (LE knot dropped).
    """
    n_total = n_thickness_cp + 1 + (1 if closed_te else 0)
    return 0.5 * (1.0 - np.cos(np.pi * np.arange(n_total) / (n_total - 1)))


def _x_camber_knots(n_camber_cp=N_CAMBER_CP):
    """Cosine-clustered x positions of the camber-spline knots."""
    return 0.5 * (1.0 - np.cos(np.pi * np.arange(n_camber_cp) / (n_camber_cp - 1)))


def _naca_half_thickness(x, tmax=0.12):
    """NACA 4-digit finite-TE half-thickness law (matches scirocco::NacaT)."""
    return (tmax / 0.2) * (
        0.2969 * np.sqrt(x)
        - 0.1260 * x
        - 0.3516 * x ** 2
        + 0.2843 * x ** 3
        - 0.1015 * x ** 4
    )


def _initial_cps(tmax=0.12):
    """NACA-00tt-equivalent initial guess: zero camber + NACA thickness law."""
    camber_cp = np.zeros(N_CAMBER_CP)
    xi_full = _xi_thickness_knots()
    # Free knots are indices 1..n_total-2 (closed TE) or 1..n_total-1 (open).
    free_slice = slice(1, 1 + N_THICKNESS_CP)
    xi_free = xi_full[free_slice]
    thickness_cp = _naca_half_thickness(xi_free ** 2, tmax=tmax)
    return camber_cp, thickness_cp


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
    """Launch both the analysis and geometry Philote servers."""
    analysis = _start_server(ANALYSIS_BINARY, ANALYSIS_ADDR)
    geometry = _start_server(GEOMETRY_BINARY, GEOMETRY_ADDR)
    # Extra grace period so both gRPC sockets are ready.
    time.sleep(1.0)
    return analysis, geometry


# ---------------------------------------------------------------------------
# Problem construction
# ---------------------------------------------------------------------------

def _build_problem(analysis_channel, geometry_channel, hdf5_file=None):
    """Wire GeometryDiscipline → SiroccoDiscipline + area through OpenMDAO."""
    nn = NPANEL + 1

    geometry_kwargs = dict(
        channel=geometry_channel,
        n_camber_cp=N_CAMBER_CP,
        n_thickness_cp=N_THICKNESS_CP,
        npanel=NPANEL,
        closed_te=CLOSED_TE,
    )
    analysis_kwargs = dict(
        channel=analysis_channel,
        naca="0012",      # initial topology used for paneling
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
        promotes_inputs=["camber_cp", "thickness_cp"],
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


def _set_baseline_inputs(prob, camber_cp, thickness_cp, alpha=5.0):
    prob.set_val("camber_cp", camber_cp)
    prob.set_val("thickness_cp", thickness_cp)
    prob.set_val("analysis.alpha", alpha)
    prob.set_val("analysis.Re", 1e6)
    prob.set_val("analysis.Ma", 0.0)


# ---------------------------------------------------------------------------
# Examples
# ---------------------------------------------------------------------------

def example_single_point(analysis_channel, geometry_channel):
    """Evaluate the spline geometry tuned to mimic NACA-0012 at alpha = 5°."""
    prob = _build_problem(analysis_channel, geometry_channel)
    prob.setup()

    camber_cp, thickness_cp = _initial_cps(tmax=0.12)
    _set_baseline_inputs(prob, camber_cp, thickness_cp, alpha=5.0)
    prob.run_model()

    cl = prob.get_val("analysis.cl")[0]
    cd = prob.get_val("analysis.cd")[0]
    cm = prob.get_val("analysis.cm")[0]
    area = prob.get_val("area_comp.area")[0]

    print("=" * 60)
    print("Single-Point Analysis: NACA-0012-equivalent spline shape")
    print("=" * 60)
    print(f"  alpha = 5.0 deg, Re = 1e6, Ma = 0.0")
    print(f"  Cl   = {cl:.6f}")
    print(f"  Cd   = {cd:.6f}")
    print(f"  Cm   = {cm:.6f}")
    print(f"  Area = {area:.6f}")
    print("=" * 60)
    return prob


def _reference_area(analysis_channel, geometry_channel, camber_cp,
                    thickness_cp):
    """Compute area for a given CP set."""
    prob = _build_problem(analysis_channel, geometry_channel)
    prob.setup()
    _set_baseline_inputs(prob, camber_cp, thickness_cp, alpha=0.0)
    prob.run_model()
    return (
        float(prob.get_val("area_comp.area")[0]),
        prob.get_val("x").copy(),
    )


def example_optimization(analysis_channel, geometry_channel, target_cl=0.5):
    """Minimize Cd subject to Cl = target_cl and area ≥ 80% of NACA-0012."""
    nn = NPANEL + 1
    hdf5_path = Path(__file__).with_name("spline_optimization_history.h5")

    camber_init, thickness_init = _initial_cps(tmax=0.12)
    init_area, x_init = _reference_area(
        analysis_channel, geometry_channel, camber_init, thickness_init
    )

    prob = _build_problem(analysis_channel, geometry_channel,
                          hdf5_file=hdf5_path)

    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options["optimizer"] = "SLSQP"
    prob.driver.options["tol"] = 1e-6
    prob.driver.options["maxiter"] = 200
    prob.driver.options["disp"] = True

    # Camber CPs are bounded to keep the airfoil sensible; thickness CPs
    # must stay positive (and the LE knot is already pinned to zero
    # internally by the discipline).
    prob.model.add_design_var(
        "camber_cp", lower=-0.10, upper=0.10, ref=0.10,
    )
    prob.model.add_design_var(
        "thickness_cp", lower=1e-3, upper=0.20, ref=0.10,
    )
    prob.model.add_design_var(
        "analysis.alpha", lower=-5.0, upper=15.0, ref=10.0,
    )

    prob.model.add_objective("analysis.cd", ref=0.01)
    prob.model.add_constraint("analysis.cl", equals=target_cl, ref=target_cl)
    min_area = 0.8 * init_area
    prob.model.add_constraint(
        "area_comp.area", lower=min_area, ref=init_area,
    )

    prob.setup()
    _set_baseline_inputs(prob, camber_init, thickness_init, alpha=5.0)

    print("\n" + "=" * 60)
    print(f"Optimization: minimize Cd, target Cl = {target_cl}")
    print(f"  DVs: {N_CAMBER_CP} camber CPs + {N_THICKNESS_CP} thickness CPs"
          f" + alpha")
    print(f"  Area constraint: area >= {min_area:.6f} "
          f"(80% of NACA-0012 init = {init_area:.6f})")
    print("=" * 60 + "\n")

    prob.run_driver()

    camber_opt = prob.get_val("camber_cp").copy()
    thickness_opt = prob.get_val("thickness_cp").copy()
    alpha = prob.get_val("analysis.alpha")[0]
    cl = prob.get_val("analysis.cl")[0]
    cd = prob.get_val("analysis.cd")[0]
    cm = prob.get_val("analysis.cm")[0]
    area = prob.get_val("area_comp.area")[0]
    x_opt = prob.get_val("x").copy()

    print("\n" + "=" * 60)
    print("Optimization Results")
    print("=" * 60)
    print(f"  camber_cp    = {np.array2string(camber_opt, precision=4)}")
    print(f"  thickness_cp = {np.array2string(thickness_opt, precision=4)}")
    print(f"  alpha = {alpha:.4f} deg")
    print(f"  Cl    = {cl:.6f}")
    print(f"  Cd    = {cd:.6f}")
    print(f"  Cm    = {cm:.6f}")
    print(f"  Area  = {area:.6f} (min = {min_area:.6f})")
    print(f"  Iteration history: {hdf5_path}")
    print("=" * 60)

    plot_results(x_init, x_opt, nn, target_cl, cd, alpha,
                 camber_opt, thickness_opt)
    return prob


def plot_results(x_init, x_opt, nn, target_cl, opt_cd, opt_alpha,
                 camber_cp, thickness_cp):
    """Plot initial vs optimized airfoils plus the spline DV values."""
    import matplotlib.pyplot as plt

    xi_init, zi_init = x_init[:nn], x_init[nn:]
    xo, zo = x_opt[:nn], x_opt[nn:]

    fig, (ax_foil, ax_cps) = plt.subplots(
        2, 1, figsize=(10, 6), gridspec_kw={"height_ratios": [3, 2]},
    )

    ax_foil.plot(xi_init, zi_init, "b--", lw=1.2, label="Initial (NACA-0012)")
    ax_foil.plot(xo, zo, "r-", lw=1.8, label="Optimized (spline)")
    ax_foil.fill_between(xo, zo, alpha=0.08, color="red")
    ax_foil.set_xlabel("x/c")
    ax_foil.set_ylabel("z/c")
    ax_foil.set_title(
        rf"Spline-Parameterized Optimization — min $C_d$ "
        rf"s.t. $C_l = {target_cl}$"
        "\n"
        rf"Optimized: $C_d = {opt_cd:.5f}$,  "
        rf"$\alpha = {opt_alpha:.2f}°$"
    )
    ax_foil.set_aspect("equal")
    ax_foil.legend(loc="upper right")
    ax_foil.grid(True, alpha=0.3)

    x_cam = _x_camber_knots()
    xi_full = _xi_thickness_knots()
    xi_free = xi_full[1:1 + N_THICKNESS_CP]
    x_thk = xi_free ** 2

    ax_cps.plot(x_cam, camber_cp, "o-", color="C0", label="camber CPs")
    ax_cps.plot(x_thk, thickness_cp, "s-", color="C3",
                label="half-thickness CPs")
    ax_cps.axhline(0.0, color="k", lw=0.5, alpha=0.5)
    ax_cps.set_xlabel("x/c (camber)  /  ξ² (thickness)")
    ax_cps.set_ylabel("CP value")
    ax_cps.legend(loc="upper right")
    ax_cps.grid(True, alpha=0.3)

    fig.tight_layout()
    out = Path(__file__).with_name("spline_optimization_result.png")
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
