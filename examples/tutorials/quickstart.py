"""Run one inviscid Scirocco analysis through the two Philote servers.

Start ``scirocco_kulfan_server`` on port 50052 and ``scirocco_server`` on
port 50051 before running this file from the repository root.
"""

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
    problem.set_val("analysis.alpha", 4.0, units="deg")
    problem.set_val("analysis.Re", 0.0)
    problem.set_val("analysis.Ma", 0.0)
    problem.run_model()

    cl = problem.get_val("analysis.cl")[0]
    cd = problem.get_val("analysis.cd")[0]
    cm = problem.get_val("analysis.cm")[0]
    print(f"Cl = {cl:.8f}")
    print(f"Cd = {cd:.3e}")
    print(f"Cm = {cm:.8f}")

    if not (0.47 < cl < 0.50 and abs(cd) < 1.0e-12):
        raise RuntimeError(
            "Unexpected quickstart result; verify the documented Philote "
            "versions and that both servers were rebuilt."
        )


if __name__ == "__main__":
    main()
