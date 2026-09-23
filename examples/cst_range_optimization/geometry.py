"""Continuous camber/thickness CST geometry for a range-based airfoil benchmark."""
from math import comb
import numpy as np
import openmdao.api as om
from scipy.special import beta


def bernstein(x, degree):
    x = np.asarray(x)
    return np.array(
        [comb(degree, i) * x**i * (1 - x) ** (degree - i) for i in range(degree + 1)]
    ).T


def stations(panels):
    if panels < 4 or panels % 2:
        raise ValueError("An even panel count >= 4 is required")
    f = np.linspace(0, 1, panels // 2 + 1)
    return np.clip(1 - 2.5 * f * (1 - f) ** 1.5 - (1 - f) ** 2.5, 0, 1)


def coordinates(thickness, camber, panels, gap):
    x = stations(panels)
    h = (
        np.sqrt(x) * (1 - x) * (bernstein(x, len(thickness) - 1) @ thickness)
        + gap * x / 2
    )
    m = x * (1 - x) * (bernstein(x, len(camber) - 1) @ camber)
    return np.r_[np.r_[x[::-1], x[1:]], np.r_[(m - h)[::-1], (m + h)[1:]]]


def area_weights(degree):
    return np.array(
        [2 * comb(degree, i) * beta(i + 1.5, degree - i + 2) for i in range(degree + 1)]
    )


def elevate(coefficients, degree):
    a = np.array(coefficients, dtype=float)
    if degree < len(a) - 1:
        raise ValueError("Degree elevation cannot reduce degree")
    while len(a) - 1 < degree:
        n = len(a)
        b = np.empty(n + 1)
        b[0], b[-1] = a[0], a[-1]
        for i in range(1, n):
            b[i] = i / n * a[i - 1] + (1 - i / n) * a[i]
        a = b
    return a


class CamberThicknessCST(om.ExplicitComponent):
    def initialize(self):
        self.options.declare("panels", default=200, types=int)
        self.options.declare("thickness_degree", default=4, types=int)
        self.options.declare("camber_degree", default=3, types=int)
        self.options.declare("gap", default=0.00126)
        self.options.declare(
            "thickness_stations", default=np.r_[np.linspace(0.15, 0.65, 21), 0.9]
        )

    def setup(self):
        n, p, q = (
            self.options[k] for k in ("panels", "thickness_degree", "camber_degree")
        )
        self.add_input("thickness_coefficients", val=np.full(p + 1, 0.17))
        self.add_input("camber_coefficients", val=np.zeros(q + 1))
        self.add_output("x", shape=2 * (n + 1), units="m")
        self.add_output("area", units="m**2")
        self.add_output(
            "thickness", shape=len(self.options["thickness_stations"]), units="m"
        )
        self.add_output("wedge", units="deg")
        x = stations(n)
        bt = (np.sqrt(x) * (1 - x))[:, None] * bernstein(x, p)
        bm = (x * (1 - x))[:, None] * bernstein(x, q)
        self.jt = np.vstack((np.zeros((n + 1, p + 1)), -bt[::-1], bt[1:]))
        self.jm = np.vstack((np.zeros((n + 1, q + 1)), bm[::-1], bm[1:]))
        self.aw = area_weights(p)
        z = self.options["thickness_stations"]
        self.tw = (2 * np.sqrt(z) * (1 - z))[:, None] * bernstein(z, p)
        self.declare_partials("x", "thickness_coefficients", val=self.jt)
        self.declare_partials("x", "camber_coefficients", val=self.jm)
        self.declare_partials(
            "area", "thickness_coefficients", val=self.aw.reshape(1, -1)
        )
        self.declare_partials("thickness", "thickness_coefficients", val=self.tw)
        for key in ("thickness_coefficients", "camber_coefficients"):
            self.declare_partials("wedge", key)

    def compute(self, inputs, outputs):
        t, m = inputs["thickness_coefficients"], inputs["camber_coefficients"]
        gap = self.options["gap"]
        outputs["x"] = coordinates(t, m, self.options["panels"], gap)
        outputs["area"] = self.aw @ t + gap / 2
        outputs["thickness"] = self.tw @ t + gap * self.options["thickness_stations"]
        upper_slope, lower_slope = -m[-1] - t[-1] + gap / 2, -m[-1] + t[-1] - gap / 2
        outputs["wedge"] = (
            (np.arctan(lower_slope) - np.arctan(upper_slope)) * 180 / np.pi
        )

    def compute_partials(self, inputs, partials):
        t, m = inputs["thickness_coefficients"], inputs["camber_coefficients"]
        gap = self.options["gap"]
        up, lo = -m[-1] - t[-1] + gap / 2, -m[-1] + t[-1] - gap / 2
        jt, jm = np.zeros(len(t)), np.zeros(len(m))
        jt[-1] = (1 / (1 + lo * lo) + 1 / (1 + up * up)) * 180 / np.pi
        jm[-1] = (-1 / (1 + lo * lo) + 1 / (1 + up * up)) * 180 / np.pi
        partials["wedge", "thickness_coefficients"] = jt
        partials["wedge", "camber_coefficients"] = jm
