"""Evaluate thickness and camber changes through the Kulfan discipline."""

import grpc
import openmdao.api as om
from philote_mdo.openmdao import RemoteExplicitComponent
from tutorial_airfoil import CST_0012, TE_UPPER, TE_LOWER


NPANEL = 100
BASE = CST_0012


def main():
    problem = om.Problem()
    problem.model.add_subsystem(
        "geometry",
        RemoteExplicitComponent(
            channel=grpc.insecure_channel("localhost:50052"),
            n_cst_upper=6,
            n_cst_lower=6,
            te_upper=TE_UPPER,
            te_lower=TE_LOWER,
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
    problem.set_val("analysis.alpha", 2.0, units="deg")
    problem.set_val("analysis.Re", 1.0e6)
    problem.set_val("analysis.Ma", 0.0)

    cases = {
        "thin": (0.80 * BASE, 0.80 * BASE),
        "baseline": (BASE, BASE),
        "cambered": (1.12 * BASE, 0.88 * BASE),
    }
    for name, (upper, lower) in cases.items():
        problem.set_val("cst_upper", upper)
        problem.set_val("cst_lower", lower)
        problem.run_model()
        print(
            f"{name:9s}",
            f"Cl={problem.get_val('analysis.cl')[0]:.5f}",
            f"Cd={problem.get_val('analysis.cd')[0]:.5f}",
        )


if __name__ == "__main__":
    main()
