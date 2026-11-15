#!/usr/bin/env python3
"""Independent 60-digit log-kernel integrals and coordinate derivatives.

Requires mpmath only when regenerating this small text fixture. Neither mfoil
nor C++ is imported: integrate L/(2*pi)*log(|xi-X(t)|) times (1-t,t), and
numerically differentiate that high-precision integral with mp.diff.
"""
from pathlib import Path
import mpmath as mp

mp.mp.dps = 60
CASES = [
    (0.0, 0.0, 0.003, 0.0002, 2.0, 0.02),
    (0.0, 0.0, 0.003, 0.0002, -2.0, -0.02),
    (0.0, 0.0, 1e-6, 2e-7, 1.0, 0.2),
    (0.0, 0.0, 0.003, 0.0002, 1000.0, 200.0),
    (2.0, -0.4, 2.003, -0.3998, 4.0, -0.38),
    (0.0, 0.0, 1.0, 0.0, 0.5, 4.001),
    (0.0, 0.0, 1.0, 0.0, 0.5, 3.999),
]


def integral(v, end):
    ax, ay, bx, by, px, py = v
    dx, dy = bx - ax, by - ay
    length = mp.sqrt(dx * dx + dy * dy)
    return (
        length
        / (2 * mp.pi)
        * mp.quad(
            lambda t: mp.log(mp.sqrt((px - ax - t * dx) ** 2 + (py - ay - t * dy) ** 2))
            * (t if end else 1 - t),
            [0, 1],
        )
    )


if __name__ == "__main__":
    output = Path(__file__).with_name("vortex_stream_reference.txt")
    with output.open("w") as f:
        f.write(
            "# ax ay bx by px py a b da/d(ax ay bx by px py) db/d(ax ay bx by px py)\n"
        )
        for case in CASES:
            v = tuple(mp.mpf(x) for x in case)
            values = [integral(v, end) for end in (0, 1)]
            deriv = []
            for end in (0, 1):
                for k in range(6):
                    order = tuple(int(j == k) for j in range(6))
                    deriv.append(mp.diff(lambda *args: integral(args, end), v, order))
            f.write(
                " ".join(format(float(x), ".17g") for x in (*v, *values, *deriv)) + "\n"
            )
    print(output)
