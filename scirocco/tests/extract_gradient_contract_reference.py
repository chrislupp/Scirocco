#!/usr/bin/env python3
"""Generate gradient references from Python mfoil, without reading C++ data.

The real oracle is bundled mfoil.py without complex-step patches. Its local
corrections are documented in the correctness review. Panel-coordinate inputs are
reparameterized using its spline2d. The CS path calls the SAME spline2d;
only SciPy's rejection of complex knots and nonanalytic primitives are
extended by mfoil_cs. Check convergence and real primal agreement before
accepting each derivative. Run test_complex_step.py before this generator.
"""
import hashlib
import importlib.util
from pathlib import Path
import sys
import numpy as np
import h5py

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'mfoil'))
import mfoil
import mfoil_cs

spec = importlib.util.spec_from_file_location('original_mfoil', ROOT/'mfoil/mfoil.py')
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)
STEP = 1e-20
OUTPUT = Path(__file__).with_name('gradient_contract_reference.h5')


def geometry(module, m):
    pp = module.spline2d(m.foil.x)
    m.foil.s = pp['X'].x.copy()
    m.foil.t = module.splinetan(pp, m.foil.s)


def model(module, naca, panels):
    m = module.mfoil(naca=naca, npanel=panels)
    m.param.verb = 0
    m.param.doplot = False
    return m


def outputs(m):
    assert m.glob.conv, 'Reference primal did not converge'
    value = np.asarray([m.post.cl, m.post.cd, m.post.cm])
    assert np.isfinite(value).all(), 'Nonfinite reference output'
    return value


def case(file, name, naca, panels, alpha, re, ma, nodes):
    base = model(original, naca, panels)
    geometry(original, base)
    kwargs = dict(alpha=alpha, Ma=ma)
    if re > 0:
        kwargs['Re'] = re
    base.setoper(**kwargs)
    base.solve()
    primal = outputs(base).real
    print(name, 'primal', primal, flush=True)
    group = file.create_group(name)
    for key, value in dict(naca=naca, panels=panels, alpha=alpha, Re=re, Ma=ma).items():
        group.attrs[key] = value
    group['x'] = base.foil.x
    group['s'] = base.foil.s
    group['outputs'] = primal
    directions = []
    for node in nodes:
        for dim in range(2):
            direction = np.zeros_like(base.foil.x)
            direction[dim, node] = 1
            directions.append(direction)
    # Global shape perturbations exercise every column, without changing topology.
    for phase in (.2, .8):
        direction = np.zeros_like(base.foil.x)
        direction[1] = .01*np.sin(np.arange(base.foil.N)*.13 + phase)
        direction[0] = .005*np.cos(np.arange(base.foil.N)*.07 + phase)
        directions.append(direction)
    jac, ds = [], []
    for i, direction in enumerate(directions):
        m = model(mfoil, naca, panels)
        m.foil.x = base.foil.x.astype(complex) + 1j*STEP*direction
        with mfoil_cs.complex_step_mode():
            geometry(mfoil, m)
            m.setoper(**kwargs)
            m.solve()
        value = outputs(m)
        np.testing.assert_allclose(value.real, primal, rtol=2e-9, atol=2e-11)
        np.testing.assert_allclose(m.foil.s.real, base.foil.s, rtol=2e-13, atol=2e-13)
        jac.append(value.imag/STEP)
        ds.append(m.foil.s.imag/STEP)
        print(name, 'shape', i, jac[-1], flush=True)
    group['directions'] = directions
    group['shape_derivatives'] = jac
    group['arc_length_derivatives'] = ds
    scalars = []
    for variable in ('alpha', 'Ma', 'Re'):
        if (variable == 'Ma' and ma == 0) or (variable == 'Re' and re == 0):
            scalars.append(np.zeros(3))
            continue
        m = model(mfoil, naca, panels)
        geometry(mfoil, m)
        perturbed = kwargs.copy()
        perturbed[variable] += 1j*STEP
        with mfoil_cs.complex_step_mode():
            m.setoper(**perturbed)
            m.solve()
        value = outputs(m)
        np.testing.assert_allclose(value.real, primal, rtol=2e-9, atol=2e-11)
        scalars.append(value.imag/STEP)
    group['scalar_derivatives'] = scalars
    file.flush()


if __name__ == '__main__':
    temp = OUTPUT.with_suffix('.tmp.h5')
    with h5py.File(temp, 'w') as f:
        for name in ('mfoil.py', 'mfoil_cs.py', 'complex_spline.py'):
            f.attrs[name+'_sha256'] = hashlib.sha256((ROOT/'mfoil'/name).read_bytes()).hexdigest()
        f.attrs['reference'] = 'Python mfoil primal; original spline2d with complex knot extension'
        case(f, 'inviscid_te', '2412', 80, 5., 0., .3, [0, 1, 40, 79, 80])
        case(f, 'viscous_spline', '0012', 100, 3., 1e6, 0., [0, 1, 20, 47, 48, 50, 80, 99, 100])
        case(f, 'viscous_mach', '0012', 100, 3., 1e6, .3, [0, 48, 50, 100])
    temp.replace(OUTPUT)
    print('Wrote', OUTPUT)
