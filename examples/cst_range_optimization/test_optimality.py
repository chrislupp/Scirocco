"""Check multiplier signs and reject feasible but nonstationary points."""
import unittest
import numpy as np
from optimality import check_kkt


class Quadratic:
    c_lower = np.array([1.0, 0.0])
    c_upper = np.array([np.inf, 0.0])
    x_lower = np.full(2, -np.inf)
    x_upper = np.full(2, np.inf)

    def _compute_objective_gradient(self, x):
        return x

    def _compute_constraints(self, x):
        return x

    def _compute_constraint_jacobian(self, x):
        return np.eye(2)


class OptimalityTest(unittest.TestCase):
    def test_lower_bound_and_equality(self):
        result = check_kkt(Quadratic(), [1.0, 0.0])
        self.assertTrue(result["passed"])
        self.assertAlmostEqual(result["multipliers"][0], 1.0)
        self.assertFalse(check_kkt(Quadratic(), [2.0, 0.0])["passed"])
        self.assertFalse(check_kkt(Quadratic(), [0.0, 0.0])["passed"])

    def test_upper_bound(self):
        model = Quadratic()
        model.c_lower = np.array([-np.inf, 0.0])
        model.c_upper = np.array([-1.0, 0.0])
        result = check_kkt(model, [-1.0, 0.0])
        self.assertTrue(result["passed"])
        self.assertAlmostEqual(result["multipliers"][0], 1.0)

    def test_nonfinite(self):
        with self.assertRaises(ValueError):
            check_kkt(Quadratic(), [np.nan, 0.0])


if __name__ == "__main__":
    unittest.main()
