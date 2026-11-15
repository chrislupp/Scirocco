"""Independent first-order KKT residuals in the optimizer's scaled coordinates."""
import numpy as np
from scipy.optimize import lsq_linear


def check_kkt(model, x, tolerance=1e-6):
    x = np.asarray(x)
    gradient = np.asarray(model._compute_objective_gradient(x))
    values = np.r_[model._compute_constraints(x), x]
    jacobian = np.vstack([model._compute_constraint_jacobian(x), np.eye(len(x))])
    if not all(np.all(np.isfinite(a)) for a in (x, gradient, values, jacobian)):
        raise ValueError("KKT check requires finite values and derivatives")
    lower = np.r_[model.c_lower, model.x_lower]
    upper = np.r_[model.c_upper, model.x_upper]
    primal = float(max(0.0, np.max(lower - values), np.max(values - upper)))
    columns = []
    multiplier_lower = []
    slacks = []
    active = []
    for i, (value, lo, hi) in enumerate(zip(values, lower, upper)):
        if lo == hi:
            columns.append(jacobian[i])
            multiplier_lower.append(-np.inf)
            slacks.append(0.0)
            active.append([i, "equality"])
        else:
            if np.isfinite(lo) and value - lo <= tolerance:
                columns.append(-jacobian[i])
                multiplier_lower.append(0.0)
                slacks.append(value - lo)
                active.append([i, "lower"])
            if np.isfinite(hi) and hi - value <= tolerance:
                columns.append(jacobian[i])
                multiplier_lower.append(0.0)
                slacks.append(hi - value)
                active.append([i, "upper"])
    if columns:
        matrix = np.column_stack(columns)
        solution = lsq_linear(
            matrix,
            -gradient,
            bounds=(multiplier_lower, np.full(len(columns), np.inf)),
            tol=1e-13,
            max_iter=1000,
        )
        multipliers = solution.x
        stationarity = float(np.max(abs(gradient + matrix @ multipliers)))
        complementarity = float(np.max(abs(multipliers * np.array(slacks))))
    else:
        multipliers = np.array([])
        stationarity = float(np.max(abs(gradient)))
        complementarity = 0.0
    return {
        "passed": max(primal, stationarity, complementarity) <= tolerance,
        "tolerance": tolerance,
        "primal_infinity_norm": primal,
        "stationarity_infinity_norm": stationarity,
        "complementarity_infinity_norm": complementarity,
        "active_constraints_and_bounds": active,
        "multipliers": multipliers.tolist(),
        "description": "Nonnegative inequality multipliers fitted to stationarity; inactive multipliers zero. All residuals use the original driver scaling.",
    }
