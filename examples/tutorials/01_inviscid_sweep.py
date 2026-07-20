"""Sweep angle of attack with the inviscid Scirocco model."""

import grpc
import numpy as np
import openmdao.api as om
from philote_mdo.openmdao import RemoteExplicitComponent


NPANEL = 100
CST_0012 = np.array([0.170, 0.160, 0.145, 0.145, 0.160, 0.170])


def main():
    problem = om.Problem()
    problem.model.add_subsystem(
        "geometry",
        RemoteExplicitComponent(
            channel=grpc.insecure_channel("localhost:50052"),
            n_cst_upper=6,
            n_cst_lower=6,
            npanel=NPANEL,
        ),
        promotes_inputs=["cst_upper", "cst_lower"],
        promotes_outputs=["x"],
    )
    problem.model.add_subsystem(
        "analysis",
        RemoteExplicitComponent(
            channel=grpc.insecure_channel("localhost:50051"),
            naca="0012",
            npanel=NPANEL,
        ),
        promotes_inputs=["x"],
    )
    problem.setup()

    problem.set_val("cst_upper", CST_0012)
    problem.set_val("cst_lower", CST_0012)
    problem.set_val("analysis.Re", 0.0)
    problem.set_val("analysis.Ma", 0.0)

    print(" alpha       Cl          Cd          Cm")
    for alpha in np.arange(-4.0, 12.1, 2.0):
        problem.set_val("analysis.alpha", alpha, units="deg")
        problem.run_model()
        cl = problem.get_val("analysis.cl")[0]
        cd = problem.get_val("analysis.cd")[0]
        cm = problem.get_val("analysis.cm")[0]
        print(f"{alpha:6.1f}  {cl:9.5f}  {cd:10.3e}  {cm:9.5f}")


if __name__ == "__main__":
    main()
