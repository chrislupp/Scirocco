"""Complex-step extension of the not-a-knot spline used by Python mfoil.

SciPy accepts complex ordinates, but rejects complex knot positions. Solve for
Hermite slopes using continuity of the second derivative and not-a-knot end
conditions instead. This is the same interpolation problem as CubicSpline;
it does not depend on Scirocco. Real inputs still use SciPy directly.
"""
import numpy as np
from scipy.interpolate import CubicSpline


class ComplexCubicSpline:
    """Scalar, not-a-knot cubic spline with complex knot positions."""

    def __init__(self, x, y):
        self.x = np.asarray(x)
        y = np.asarray(y)
        n = len(x)
        if n < 2 or np.any(np.diff(self.x.real) <= 0):
            raise ValueError("Spline knots must be strictly increasing")
        h = np.diff(self.x)
        secant = np.diff(y) / h
        dtype = np.result_type(x, y, complex)
        if n == 2:
            slopes = np.full(2, secant[0], dtype=dtype)
        else:
            a = np.zeros((n, n), dtype=dtype)
            b = np.zeros(n, dtype=dtype)
            for i in range(1, n - 1):
                a[i, i-1:i+2] = h[i], 2*(h[i-1]+h[i]), h[i-1]
                b[i] = 3*(h[i]*secant[i-1] + h[i-1]*secant[i])
            if n == 3:
                # Both not-a-knot conditions reduce to a single parabola.
                a[0, :2] = a[-1, -2:] = (1, 1)
                b[0], b[-1] = 2*secant[0], 2*secant[-1]
            else:
                d = h[0] + h[1]
                a[0, :2] = h[1], d
                b[0] = ((h[0]+2*d)*h[1]*secant[0] + h[0]**2*secant[1])/d
                d = h[-2] + h[-1]
                a[-1, -2:] = d, h[-2]
                b[-1] = (h[-1]**2*secant[-2] + (2*d+h[-1])*h[-2]*secant[-1])/d
            slopes = np.linalg.solve(a, b)
        # Hermite interpolation on each interval, in descending powers.
        t = (slopes[:-1] + slopes[1:] - 2*secant)/h
        self.c = np.array([t/h, (secant-slopes[:-1])/h-t, slopes[:-1], y[:-1]])

    def __call__(self, values):
        values = np.asarray(values)
        i = np.clip(np.searchsorted(self.x.real, values.real, side='right')-1,
                    0, len(self.x)-2)
        t = values-self.x[i]
        result = np.zeros_like(t, dtype=np.result_type(self.c, t))
        for c in self.c:
            result = result*t + c[i]
        return result

    def derivative(self):
        result = object.__new__(ComplexCubicSpline)
        result.x = self.x
        degree = len(self.c)-1
        result.c = self.c[:-1]*np.arange(degree, 0, -1)[:, None]
        return result


def cs_cubic_spline(x, y):
    if np.iscomplexobj(x):
        return ComplexCubicSpline(x, y)
    return CubicSpline(x, y)
