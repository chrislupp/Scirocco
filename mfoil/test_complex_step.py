"""Validate the CS extension against Python mfoil and interpolation identities.

Run: python3 -m unittest discover -s mfoil -p 'test_*.py'
No C++ values or generated C++ data are used here.
"""
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
