"""Write a viscous angle-of-attack sweep to a server-side HDF5 history."""

import argparse

import grpc
import numpy as np
import openmdao.api as om
from philote_mdo.openmdao import RemoteExplicitComponent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--history-file", default="analysis_history.h5")
    args = parser.parse_args()

    npanel = 100
    base = np.array([0.170, 0.160, 0.145, 0.145, 0.160, 0.170])
    problem = om.Problem()
    problem.model.add_subsystem(
        "geometry",
        RemoteExplicitComponent(
            channel=grpc.insecure_channel("localhost:50052"),
            n_cst_upper=6,
            n_cst_lower=6,
            npanel=npanel,
        ),
        promotes_inputs=["cst_upper", "cst_lower"],
        promotes_outputs=["x"],
    )
    problem.model.add_subsystem(
        "analysis",
        RemoteExplicitComponent(
            channel=grpc.insecure_channel("localhost:50051"),
            naca="0012",
            npanel=npanel,
            hdf5_file=args.history_file,
        ),
        promotes_inputs=["x"],
    )
    problem.setup()
    problem.set_val("cst_upper", base)
    problem.set_val("cst_lower", base)
    problem.set_val("analysis.Re", 1.0e6)
    problem.set_val("analysis.Ma", 0.0)

    for alpha in np.linspace(-2.0, 8.0, 11):
        problem.set_val("analysis.alpha", alpha, units="deg")
        problem.run_model()
        print(
            f"alpha={alpha:5.1f}",
            f"Cl={problem.get_val('analysis.cl')[0]:.5f}",
            f"Cd={problem.get_val('analysis.cd')[0]:.5f}",
        )

    print(
        f"Stop scirocco_server to flush {args.history_file}, then run "
        "examples/tutorials/read_hdf5_history.py."
    )


if __name__ == "__main__":
    main()
