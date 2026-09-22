"""Validate the CS extension against Python mfoil and interpolation identities.

Run: python3 -m unittest discover -s mfoil -p 'test_*.py'
No C++ values or generated C++ data are used here.
"""
import copy
import importlib.util
from pathlib import Path
import unittest
import numpy as np
from numpy.testing import assert_allclose
from scipy.interpolate import CubicSpline
from complex_spline import ComplexCubicSpline
import mfoil_cs
import mfoil

spec = importlib.util.spec_from_file_location('original_mfoil', Path(__file__).with_name('mfoil.py'))
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)


class ComplexStepTest(unittest.TestCase):
    def test_pressure_coefficients_and_partials_are_scale_invariant(self):
        base = original.mfoil(naca='2412', npanel=40)
        base.param.verb = 0
        base.setoper(alpha=4.0, Ma=0.3)
        base.solve()
        expected = np.array([base.post.cl, base.post.cm, base.post.cdpi])
        # The same pressure field on geometrically similar foils gives the
        # same dimensionless forces, moments, and velocity/angle partials.
        for scale in (0.5, 2.0, 3.2):
            for module in (original, mfoil):
                model = copy.deepcopy(base)
                model.foil.x *= scale
                model.geom.xref *= scale
                model.geom.chord *= scale
                module.calc_force(model)
                assert_allclose([model.post.cl, model.post.cm, model.post.cdpi],
                                expected, rtol=2e-13, atol=2e-15)
                assert_allclose(model.post.cl_ue, base.post.cl_ue,
                                rtol=2e-13, atol=2e-15)
                assert_allclose(model.post.cl_alpha, base.post.cl_alpha,
                                rtol=2e-13, atol=2e-15)

    def test_nonunit_pressure_lift_partials_against_complex_step(self):
        base = original.mfoil(naca='2412', npanel=40)
        base.param.verb = 0
        base.foil.x *= 2.0
        base.geom.chord = 2.0
        base.setoper(alpha=4.0, Ma=0.3)
        base.solve()
        step = 1e-20
        for node in range(base.foil.N):
            model = copy.deepcopy(base)
            # Perturb one physical edge velocity, holding the others fixed.
            model.isol.gamref = model.isol.gamref.astype(complex)
            model.isol.gamref[node,0] += 1j*step/original.cosd(base.oper.alpha)
            with mfoil_cs.complex_step_mode():
                mfoil.calc_force(model)
            assert_allclose(model.post.cl.imag/step, base.post.cl_ue[node],
                            rtol=2e-13, atol=2e-15)
        model = copy.deepcopy(base)
        model.oper.alpha += 1j*step
        with mfoil_cs.complex_step_mode():
            mfoil.calc_force(model)
        due = (-original.sind(base.oper.alpha)*base.isol.gamref[:,0]
               + original.cosd(base.oper.alpha)*base.isol.gamref[:,1])*np.pi/180
        expected = base.post.cl_alpha + base.post.cl_ue @ due
        assert_allclose(model.post.cl.imag/step, expected, rtol=2e-13, atol=2e-15)

    def test_nonunit_viscous_lift_trim(self):
        models = []
        for module, target in ((original, 0.4), (mfoil, 0.4+1j*1e-20)):
            model = module.mfoil(naca='0012', npanel=100)
            model.param.verb = 0
            model.foil.x *= 2.0
            model.geom.chord = 2.0
            model.foil.s = module.spline2d(model.foil.x)['X'].x.copy()
            if module is original:
                model.setoper(alpha=3.0, Re=1e6, Ma=0.3, cl=target)
                model.solve()
            else:
                with mfoil_cs.complex_step_mode():
                    model.setoper(alpha=3.0, Re=1e6, Ma=0.3, cl=target)
                    model.solve()
            self.assertTrue(model.glob.conv)
            assert_allclose(model.post.cl.real, 0.4, rtol=0, atol=2e-12)
            models.append(model)
        assert_allclose(models[1].oper.alpha.real, models[0].oper.alpha,
                        rtol=0, atol=2e-11)
        assert_allclose(models[1].post.cl.imag/1e-20, 1.0, rtol=0, atol=2e-12)
        self.assertTrue(np.isfinite(models[1].oper.alpha.imag/1e-20))
        self.assertGreater(abs(models[1].oper.alpha.imag/1e-20), 1.0)

    def test_full_gradient_against_high_precision_python(self):
        reference = Path(__file__).resolve().parents[1]/'scirocco/tests/high_precision_gradient_reference.txt'
        model = mfoil.mfoil(naca='0012', npanel=100)
        model.param.verb = 0
        model.param.doplot = False
        model.geom.chord = 2.0
        model.foil.x = (2*model.foil.x).astype(complex)
        model.foil.x[0,50] += 1j*1e-20
        with mfoil_cs.complex_step_mode():
            pp = mfoil.spline2d(model.foil.x)
            model.foil.s = pp['X'].x.copy()
            model.foil.t = mfoil.splinetan(pp, model.foil.s)
            model.setoper(alpha=3.0, Re=1e6, Ma=0.3)
            model.solve()
            self.assertTrue(model.glob.conv)
            values = np.array([model.post.cl, model.post.cd,
                               model.post.cm, model.post.cdpi])
        assert_allclose(values.imag/1e-20, np.loadtxt(reference), rtol=0, atol=2e-13)

    def test_vortex_stream_against_high_precision_integrals(self):
        reference = Path(__file__).resolve().parents[1]/'scirocco/tests/vortex_stream_reference.txt'
        for row in np.loadtxt(reference):
            # Fixture integrates the log-distance kernel with mpmath at 60
            # digits and differentiates that integral independently.
            coords, expected = row[:6], row[6:8]
            xj = coords[:4].reshape(2, 2).T
            xi = coords[4:]
            length = np.linalg.norm(xj[:,1]-xj[:,0])
            for module in (original, mfoil):
                got = module.panel_linvortex_stream(xj, xi)
                assert_allclose(got, expected, rtol=2e-13, atol=5e-16*length)
            for k in range(6):
                perturbed = coords.astype(complex)
                perturbed[k] += 1j*1e-20
                with mfoil_cs.complex_step_mode():
                    got = mfoil.panel_linvortex_stream(
                        perturbed[:4].reshape(2, 2).T, perturbed[4:])
                derivative = np.imag(got)/1e-20
                assert_allclose(derivative, row[[8+k, 14+k]], rtol=2e-13, atol=5e-14)

    def test_source_stream_against_high_precision_integrals(self):
        reference = Path(__file__).resolve().parents[1]/'scirocco/tests/source_stream_reference.txt'
        def values(module, coords):
            xj = coords[:4].reshape(2, 2).T
            xi = coords[4:]
            return np.array([module.panel_constsource_stream(xj, xi),
                             *module.panel_linsource_stream(xj, xi)])
        for row in np.loadtxt(reference):
            coords = row[:6]
            length = np.linalg.norm(coords[2:4]-coords[:2])
            for module in (original, mfoil):
                assert_allclose(values(module, coords), row[6:9],
                                rtol=2e-13, atol=5e-16*length)
            for k in range(6):
                perturbed = coords.astype(complex)
                perturbed[k] += 1j*1e-20
                with mfoil_cs.complex_step_mode():
                    derivative = values(mfoil, perturbed).imag/1e-20
                assert_allclose(derivative, row[[9+k,15+k,21+k]],
                                rtol=2e-13, atol=5e-14)

    def test_spline_real_values_match_scipy(self):
        for n in (2, 3, 8):
            x = np.linspace(0, 1, n)**1.5
            y = np.sin(3*x)
            ref = CubicSpline(x, y)
            cs = ComplexCubicSpline(x.astype(complex), y.astype(complex))
            assert_allclose(cs.c, ref.c, rtol=3e-12, atol=3e-12)
            t = np.linspace(-.1, 1.1, 29)
            assert_allclose(cs(t), ref(t), rtol=2e-12, atol=2e-12)
            assert_allclose(cs.derivative()(t), ref.derivative()(t), atol=2e-12)

    def test_moving_knots_reproduce_polynomial(self):
        # Interpolating a cubic at moving knots must reproduce the SAME cubic.
        # This checks the knot-position chain independently of either solver.
        for n, degree in ((2, 1), (3, 2), (8, 3)):
            x = np.linspace(0, 1, n)**1.3 + 1j*1e-20*np.cos(np.arange(n))
            y = x**degree + 2*x + 1
            p = ComplexCubicSpline(x, y)
            t = np.linspace(.05, .95, 13) + 1j*1e-20
            expected = degree*t.real**(degree-1) + 2
            assert_allclose(p(t).imag/1e-20, expected, rtol=3e-12, atol=3e-12)

    def test_spline2d_preserves_original_python_primal(self):
        for naca in ('0012', '2412'):
            foil = original.mfoil(naca=naca, npanel=60)
            ref = original.spline2d(foil.foil.x)
            with mfoil_cs.complex_step_mode():
                cs = mfoil.spline2d(foil.foil.x.astype(complex))
            for dim in ('X', 'Y'):
                assert_allclose(cs[dim].x.real, ref[dim].x, atol=2e-14, rtol=2e-14)
                assert_allclose(cs[dim].c.real, ref[dim].c, atol=2e-8, rtol=2e-8)
            # Uniform scaling multiplies every arc length by the same factor.
            with mfoil_cs.complex_step_mode():
                scaled = mfoil.spline2d(foil.foil.x*(1+1j*1e-20))
            assert_allclose(scaled['X'].x.imag/1e-20, ref['X'].x, rtol=2e-12, atol=2e-12)

    def test_update_preserves_python_projections_and_wake(self):
        base = original.mfoil(naca='0012', npanel=100)
        base.param.verb = 0
        base.param.doplot = False
        base.setoper(alpha=3., Re=1e6, Ma=.3)
        base.solve()
        self.assertTrue(base.glob.conv)
        for redowake in (False, True):
            real = copy.deepcopy(base)
            real.oper.redowake = redowake
            real.glob.dU = np.zeros_like(real.glob.U)
            real.glob.dalpha = .1
            # Force a post-update Hk correction, which the old wrapper omitted.
            real.glob.U[1, 20] = .5 * real.glob.U[0, 20]
            cs = copy.deepcopy(real)
            original.update_state(real)
            mfoil_cs._patched_update_state(cs)
            assert_allclose(cs.glob.U, real.glob.U, rtol=2e-13, atol=2e-13)
            assert_allclose(cs.isol.gam, real.isol.gam, rtol=2e-13, atol=2e-13)
            assert_allclose(cs.wake.x, real.wake.x, rtol=2e-13, atol=2e-13)

    def test_trailing_edge_absolute_value_chain(self):
        x = original.mfoil(naca='2412', npanel=80).foil.x
        real = original.TE_info(x)
        assert_allclose(mfoil.TE_info(x)[3], real[3], atol=1e-15)
        # Independently differentiate normalize(v) = (I-tt^T) dv / |v|.
        direction = np.zeros_like(x); direction[0, 0] = 1
        def normalized(v, dv):
            t = v/np.linalg.norm(v)
            return t, (dv-t*np.dot(t, dv))/np.linalg.norm(v)
        t1, dt1 = normalized(x[:, 0]-x[:, 1], direction[:, 0]-direction[:, 1])
        t2, dt2 = normalized(x[:, -1]-x[:, -2], direction[:, -1]-direction[:, -2])
        t, dt = normalized(t1+t2, dt1+dt2)
        p, dp = normalized(x[:, -1]-x[:, 0], direction[:, -1]-direction[:, 0])
        cross = t[0]*p[1]-t[1]*p[0]
        expected = np.sign(cross)*(dt[0]*p[1]+t[0]*dp[1]-dt[1]*p[0]-t[1]*dp[0])
        got = mfoil.TE_info(x+1j*1e-20*direction)[3].imag/1e-20
        assert_allclose(got, expected, rtol=1e-12)
        self.assertGreater(abs(expected), .01)  # detects imaginary-part deletion


if __name__ == '__main__':
    unittest.main()
