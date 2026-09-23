#!/usr/bin/env python3
"""Independent 60-digit source-angle integrals and coordinate derivatives.

Requires mpmath only to regenerate the fixture. Evaluates the defining angle
integral and mp.diff; imports neither mfoil nor the C++ implementation.
"""
from pathlib import Path
import mpmath as mp
from extract_vortex_stream_reference import CASES

mp.mp.dps = 60
# Unlike the log kernel, the angle kernel has a branch cut along the panel
# extension. Keep the short-panel case off that cut; test both sides below.
SOURCE_CASES = (
    CASES[:2]
    + [(0.0, 0.0, 1e-6, 2e-7, 1.0, 0.23)]
    + CASES[3:5]
    + [
        (0.0, 0.0, 1.0, 0.0, 0.4, 4.001),
        (0.0, 0.0, 1.0, 0.0, 0.4, 3.99),
        (0.0, 0.0, 1.0, 0.0, 5.0, 1e-5),
        (0.0, 0.0, 1.0, 0.0, 5.0, -1e-5),
        (0.0, 0.0, 1.0, 0.0, 0.49999, 5.0),
        (0.0, 0.0, 1.0, 0.0, 0.50001, 5.0),
    ]
)


def integral(v, kind):
    ax, ay, bx, by, px, py = v
    dx, dy = bx - ax, by - ay
    length = mp.sqrt(dx * dx + dy * dy)
    x = ((px - ax) * dx + (py - ay) * dy) / length
    z = (-(px - ax) * dy + (py - ay) * dx) / length
    offset = (
        mp.mpf("-.25")
        if mp.atan2(z, x) + mp.atan2(z, x - length) > mp.pi
        else mp.mpf(".75")
    )

    def integrand(t):
        theta = mp.atan2(z, x - t * length)
        if kind == 0:
            return theta / (2 * mp.pi) + offset
        return (theta % (2 * mp.pi)) / (2 * mp.pi) * (1 - t if kind == 1 else t)

    return length * mp.quad(integrand, [0, 1])


if __name__ == "__main__":
    output = Path(__file__).with_name("source_stream_reference.txt")
    with output.open("w") as f:
        f.write(
            "# ax ay bx by px py constant a b; six coordinate derivatives for each output\n"
        )
        for case in SOURCE_CASES:
            v = tuple(mp.mpf(x) for x in case)
            values = [integral(v, kind) for kind in range(3)]
            derivatives = []
            for kind in range(3):
                for k in range(6):
                    order = tuple(int(j == k) for j in range(6))
                    derivatives.append(mp.diff(lambda *a: integral(a, kind), v, order))
            f.write(
                " ".join(format(float(x), ".17g") for x in (*v, *values, *derivatives))
                + "\n"
            )
    print(output)
