"""Compare viscous Scirocco solutions across Reynolds number."""

import grpc
import openmdao.api as om
from philote_mdo.openmdao import RemoteExplicitComponent
from tutorial_airfoil import CST_0012, TE_UPPER, TE_LOWER


NPANEL = 100



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

    problem.set_val("cst_upper", CST_0012)
    problem.set_val("cst_lower", CST_0012)
    problem.set_val("analysis.alpha", 2.0, units="deg")
    problem.set_val("analysis.Ma", 0.0)

    print("       Re       Cl        Cd        Cm")
    for reynolds in (2.5e5, 5.0e5, 1.0e6, 2.0e6):
        problem.set_val("analysis.Re", reynolds)
        problem.run_model()
        cl = problem.get_val("analysis.cl")[0]
        cd = problem.get_val("analysis.cd")[0]
        cm = problem.get_val("analysis.cm")[0]
        print(f"{reynolds:9.2e}  {cl:8.5f}  {cd:8.5f}  {cm:8.5f}")


if __name__ == "__main__":
    main()
