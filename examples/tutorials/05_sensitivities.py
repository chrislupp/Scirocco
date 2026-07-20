"""Check the analytic CST-to-aerodynamics total derivatives."""

import grpc
import numpy as np
import openmdao.api as om
from philote_mdo.openmdao import RemoteExplicitComponent


def _active_error(error):
    return error.forward if error.forward is not None else error.reverse


def main():
    npanel = 100
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
        ),
        promotes_inputs=["x"],
    )
    for name in ("cst_upper", "cst_lower", "analysis.alpha"):
        problem.model.add_design_var(name)
    problem.model.add_objective("analysis.cd")
    problem.setup()

    base = np.array([0.170, 0.160, 0.145, 0.145, 0.160, 0.170])
    problem.set_val("cst_upper", base)
    problem.set_val("cst_lower", base)
    problem.set_val("analysis.alpha", 4.0, units="deg")
    problem.set_val("analysis.Re", 0.0)
    problem.set_val("analysis.Ma", 0.0)
    problem.run_model()

    totals = problem.check_totals(
        of=["analysis.cl", "analysis.cd", "analysis.cm"],
        wrt=["cst_upper", "cst_lower", "analysis.alpha"],
        method="fd",
        step=1.0e-6,
        compact_print=True,
    )

    absolute_errors = []
    relative_errors = []
    for metadata in totals.values():
        absolute = _active_error(metadata["abs error"])
        relative = _active_error(metadata["rel error"])
        if absolute is not None and np.isfinite(absolute):
            absolute_errors.append(absolute)
        if relative is not None and np.isfinite(relative):
            relative_errors.append(relative)

    max_absolute = max(absolute_errors)
    max_relative = max(relative_errors)
    print(f"maximum absolute error: {max_absolute:.3e}")
    print(f"maximum finite relative error: {max_relative:.3e}")
    if max_absolute > 1.0e-4 or max_relative > 1.0e-4:
        raise RuntimeError("analytic total-derivative check exceeded tolerance")


if __name__ == "__main__":
    main()
