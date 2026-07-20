"""Supply normalized NACA 0012 panel nodes from Python."""

import grpc
import numpy as np
import openmdao.api as om
from philote_mdo.openmdao import RemoteExplicitComponent


def naca0012_coordinates(npanel=100, thickness=0.12):
    if npanel < 4 or npanel % 2:
        raise ValueError("npanel must be an even integer of at least 4")

    nhalf = npanel // 2
    beta = np.linspace(0.0, np.pi, nhalf + 1)
    x_lower = 0.5 * (1.0 + np.cos(beta))
    x_upper = x_lower[-2::-1]

    def half_thickness(x):
        return (thickness / 0.2) * (
            0.2969 * np.sqrt(x)
            - 0.1260 * x
            - 0.3516 * x**2
            + 0.2843 * x**3
            - 0.1015 * x**4
        )

    z_lower = -half_thickness(x_lower)
    z_upper = +half_thickness(x_upper)
    return np.concatenate([x_lower, x_upper, z_lower, z_upper])


def main():
    npanel = 100
    coordinates = naca0012_coordinates(npanel=npanel)
    if coordinates.size != 2 * (npanel + 1):
        raise RuntimeError("coordinate vector has the wrong size")

    problem = om.Problem()
    ivc = problem.model.add_subsystem(
        "inputs", om.IndepVarComp(), promotes_outputs=["x"]
    )
    ivc.add_output("x", val=coordinates, units="m")
    problem.model.add_subsystem(
        "analysis",
        RemoteExplicitComponent(
            channel=grpc.insecure_channel("localhost:50051"),
            naca="0012",
            npanel=npanel,
        ),
        promotes_inputs=["x"],
    )
    problem.setup()
    problem.set_val("analysis.alpha", 4.0, units="deg")
    problem.set_val("analysis.Re", 1.0e6)
    problem.set_val("analysis.Ma", 0.0)
    problem.run_model()

    print(f"Cl = {problem.get_val('analysis.cl')[0]:.8f}")
    print(f"Cd = {problem.get_val('analysis.cd')[0]:.8f}")
    print(f"Cm = {problem.get_val('analysis.cm')[0]:.8f}")


if __name__ == "__main__":
    main()
