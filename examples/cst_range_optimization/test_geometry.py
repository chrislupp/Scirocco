"""Geometry derivative and representation-invariance checks (no flow servers)."""
import unittest
import numpy as np
import openmdao.api as om
from openmdao.utils.assert_utils import assert_check_partials
from scipy.integrate import quad
from geometry import CamberThicknessCST, area_weights, bernstein, coordinates, elevate


class GeometryTest(unittest.TestCase):
    def test_derivatives_at_cambered_shape(self):
        p = om.Problem()
        p.model.add_subsystem("shape", CamberThicknessCST(), promotes=["*"])
        p.setup(force_alloc_complex=True)
        p.set_val("thickness_coefficients", [0.13, 0.18, 0.12, 0.2, 0.09])
        p.set_val("camber_coefficients", [0.02, 0.08, -0.01, 0.04])
        p.run_model()
        result = p.check_partials(method="cs", out_stream=None)
        assert_check_partials(result, atol=2e-12, rtol=2e-12)
        p.cleanup()

    def test_area_and_degree_elevation(self):
        t = np.array([0.13, 0.18, 0.12, 0.2, 0.09])
        m = np.array([0.02, 0.08, -0.01, 0.04])
        gap = 0.00126
        integral = quad(
            lambda x: float(
                2 * np.sqrt(x) * (1 - x) * (bernstein([x], 4) @ t)[0] + gap * x
            ),
            0,
            1,
            epsabs=1e-13,
            epsrel=1e-13,
        )[0]
        self.assertAlmostEqual(
            integral, float(area_weights(4) @ t + gap / 2), places=13
        )
        for panels in [100, 200, 400]:
            np.testing.assert_allclose(
                coordinates(t, m, panels, gap),
                coordinates(elevate(t, 6), elevate(m, 5), panels, gap),
                atol=2e-16,
                rtol=0,
            )
        self.assertAlmostEqual(
            float(area_weights(6) @ elevate(t, 6)),
            float(area_weights(4) @ t),
            places=14,
        )


if __name__ == "__main__":
    unittest.main()
