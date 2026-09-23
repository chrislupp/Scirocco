#!/usr/bin/env python3
"""Independent Python checks for fixed arc labels, trim, and pressure drag.

The bundled Python contains the historical corrections listed in the upstream
provenance audit. Its real module is loaded without complex-step monkey patches.
The inviscid trim below solves Cl(alpha)=target around that Python model; the
published Python leaves its inviscid trim call commented out.
"""
from pathlib import Path
import hashlib
import numpy as np
import h5py
import extract_gradient_contract_reference as ref

OUTPUT = Path(__file__).with_name("remaining_gradient_reference.h5")


def solve(module, m, kwargs, redowake):
    m.oper.redowake = redowake
    m.setoper(**kwargs)
    if "cl" in kwargs and not m.oper.viscous:
        # Scalar Newton on the Python pressure integral, including its analytic
        # gamma rotation. This does not use any C++ results or sensitivities.
        for _ in range(50):
            module.solve_inviscid(m)
            residual = m.post.cl - m.oper.cltgt
            g = (
                (
                    -module.sind(m.oper.alpha) * m.isol.gamref[:, 0]
                    + module.cosd(m.oper.alpha) * m.isol.gamref[:, 1]
                )
                * np.pi
                / 180
            )
            slope = m.post.cl_alpha + m.post.cl_ue @ g
            m.oper.alpha -= residual / slope
            if abs(residual.real) < 1e-13:
                module.solve_inviscid(m)
                break
        assert abs((m.post.cl - m.oper.cltgt).real) < 1e-12
    else:
        m.solve()
    assert m.glob.conv
    values = np.asarray([m.post.cl, m.post.cd, m.post.cm, m.post.cdpi, m.oper.alpha])
    assert np.isfinite(values).all()
    return values


def case(
    file,
    name,
    naca="0012",
    alpha=3.0,
    re=1e6,
    ma=0.3,
    mode="spline",
    target=None,
    redowake=False,
    chord=1.0,
):
    panels = 100

    def model(module):
        m = ref.model(module, naca, panels)
        m.foil.x *= chord
        m.foil.s *= chord
        m.geom.chord = chord
        return m

    base = model(ref.original)
    if mode == "spline":
        ref.geometry(ref.original, base)
    kwargs = dict(alpha=alpha, Ma=ma)
    if re:
        kwargs["Re"] = re
    if target is not None:
        kwargs["cl"] = target
    primal = solve(ref.original, base, kwargs, redowake).real
    print(name, primal, flush=True)
    group = file.create_group(name)
    for key, value in dict(
        naca=naca,
        alpha=alpha,
        Re=re,
        Ma=ma,
        mode=mode,
        target=-1 if target is None else target,
        redowake=int(redowake),
        chord=chord,
    ).items():
        group.attrs[key] = value
    group["x"] = base.foil.x
    group["s"] = base.foil.s
    group["outputs"] = primal
    directions = []
    for node in (0, 48, 50, 100):
        for d in range(2):
            v = np.zeros_like(base.foil.x)
            v[d, node] = 1
            directions.append(v)
    for phase in (0.2, 0.8):
        v = np.empty_like(base.foil.x)
        v[0] = 0.005 * np.cos(np.arange(base.foil.N) * 0.07 + phase)
        v[1] = 0.01 * np.sin(np.arange(base.foil.N) * 0.13 + phase)
        directions.append(v)
    derivatives = []
    for i, direction in enumerate(directions):
        m = model(ref.mfoil)
        m.foil.x = base.foil.x.astype(complex) + 1j * ref.STEP * direction
        with ref.mfoil_cs.complex_step_mode():
            if mode == "spline":
                ref.geometry(ref.mfoil, m)
            value = solve(ref.mfoil, m, kwargs, redowake)
        np.testing.assert_allclose(value.real, primal, atol=3e-10, rtol=3e-9)
        derivatives.append(value.imag / ref.STEP)
        print(name, "shape", i, derivatives[-1], flush=True)
    group["directions"] = directions
    group["shape_derivatives"] = derivatives
    scalars = []
    for variable in ("alpha", "Ma", "Re", "cl"):
        if variable not in kwargs or (variable == "Ma" and ma == 0):
            scalars.append(np.zeros(5))
            continue
        m = model(ref.mfoil)
        if mode == "spline":
            ref.geometry(ref.mfoil, m)
        kw = kwargs.copy()
        kw[variable] += 1j * ref.STEP
        with ref.mfoil_cs.complex_step_mode():
            value = solve(ref.mfoil, m, kw, redowake)
        np.testing.assert_allclose(value.real, primal, atol=3e-10, rtol=3e-9)
        scalars.append(value.imag / ref.STEP)
        print(name, variable, scalars[-1], flush=True)
    group["scalar_derivatives"] = scalars
    file.flush()


if __name__ == "__main__":
    temp = OUTPUT.with_suffix(".tmp.h5")
    with h5py.File(temp, "w") as f:
        for name in ("mfoil.py", "mfoil_cs.py", "complex_spline.py"):
            f.attrs[name + "_sha256"] = hashlib.sha256(
                (ref.ROOT / "mfoil" / name).read_bytes()
            ).hexdigest()
        case(f, "fixed_arc_labels", mode="fixed")
        case(f, "negative_alpha", alpha=-3.0)
        case(f, "inviscid_trim", naca="2412", re=0, target=0.5)
        case(f, "viscous_trim_frozen", target=0.4)
        case(f, "viscous_trim_rebuilt", target=0.4, redowake=True)
        case(f, "nonunit_chord", chord=2.0)
        case(f, "cambered_viscous", naca="2412", alpha=2.0, ma=0.2)
    temp.replace(OUTPUT)
