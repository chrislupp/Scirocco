Fixed-panel multipoint airfoil optimization
==========================================

Run the example from the repository root, using the Python environment with
OpenMDAO, modopt/IPOPT, philote_mdo, NumPy, SciPy, and grpc installed::

    python examples/cst_range_optimization/run.py \
        --binary BUILD/src/philote/scirocco_server --output OUTPUT/multipoint

The default is one optimization at 100 panels, initialized from a symmetric
NACA 0012 CST fit. It uses analytic derivatives, IPOPT with 20 stored BFGS
updates, and at most 160 iterations. It does not launch panel studies or
mesh-refinement stages. A fresh output directory is required for each run.

A successful result must satisfy both the optimizer's success criterion
and an independent first-order KKT check at tolerance 1e-6, using the same
scaled variables and constraints. A failed trial flow analysis is reported
to IPOPT for step rejection; no artificial lift, drag, or gradients replace it.
The run saves its design, forces, convergence residuals, source copies and
hashes, evaluation history, and server process manifest. Servers are closed
when the run exits.

Problem
-------

Re = 1e6, Mach = 0, and Ncrit = 9 are fixed. Six lift targets contribute
to weighted drag: Cl = 0.3, 0.4, 0.5, 0.6, 0.7, and 0.8. Two additional
conditions at Cl = 0.25 and 0.85 constrain off-design behavior. Eight angle
variables enforce the lift targets.

Nine shape variables describe continuous camber and half-thickness::

    h(x) = sqrt(x) (1-x) T(x) + gap*x/2
    m(x) = x (1-x) M(x)
    upper = m+h; lower = m-h

T has five Bernstein coefficients and M has four. Both surfaces share
leading-edge radius T(0)^2/2; the trailing-edge gap is fixed at 0.00126.
The objective weights are 0.08, 0.20, 0.24, 0.2133333333, 0.1866666667,
and 0.08. Guard conditions have zero objective weight.

The following are explicit choices for a general benchmark, not universal
airfoil design requirements:

* Leading-edge radius between 0.005 and 0.020 chord.
* Cm >= -0.08 at every operating point.
* Continuous section area >= 80% of the reference area.
* Thickness >= 80% of reference at 21 stations from 0.15 to 0.65 chord,
  and at 0.90 chord.
* Trailing-edge wedge between 4 and 20 degrees.
* Guard-condition drag <= 105% of matched reference drag.
* Angles between -2 and 10 degrees; thickness coefficients nonnegative.
* No coefficient-curvature constraint.

The reference is the unchanged continuous NACA 0012 CST fit used by the
existing multipoint example. The degree-4 initial thickness distribution is
its least-squares approximation. Both reference and initial forces are saved.

Convergence checks
------------------

Two documented, modestly perturbed starting shapes are available::

    python examples/cst_range_optimization/run.py --binary BUILD/src/philote/scirocco_server \
        --start cambered --output OUTPUT/cambered
    python examples/cst_range_optimization/run.py --binary BUILD/src/philote/scirocco_server \
        --start thicker --output OUTPUT/thicker

The cambered start adds 0.5% chord peak camber and a smooth +/-3% thickness
coefficient perturbation. The thicker start increases thickness coefficients
by 3%. Neither uses an optimized design. Success from these finite tests is
evidence of practical reliability, not a guarantee of global convergence.

Geometry and KKT-checker tests require no flow servers::

    python examples/cst_range_optimization/test_geometry.py
    python examples/cst_range_optimization/test_optimality.py

The default run checks aerodynamic directional derivatives at three central
finite-difference steps. --skip-check skips only this diagnostic, not analytic
gradients or final convergence checks. To check every shape coefficient and
angle derivative at a saved design::

    python examples/cst_range_optimization/check_gradients.py OUTPUT/multipoint/summary.json \
        --binary BUILD/src/philote/scirocco_server --output OUTPUT/gradients

Results
-------

Plot the geometry and design-point forces::

    python examples/cst_range_optimization/plot.py OUTPUT/multipoint/summary.json \
        --output OUTPUT/multipoint/result

summary.json distinguishes solver success from the independent KKT check.
A run stopped by its iteration limit is not labeled an optimum. Optional
research utilities remain available in validate.py; they are not part of the
default example. --seed permits an explicit continuation; unfinished seeds
require --allow-unconverged-seed and retain their original status in the log.

Background
----------

The camber/thickness representation and multipoint formulation are motivated
by these sources; the numerical bounds above are benchmark choices:

* https://web.mit.edu/drela/Public/papers/Pros_Cons_Airfoil_Optimization.pdf
* https://www.brendakulfan.com/_files/ugd/169bff_ccc8f39379d444298b9336f8629f49cb.pdf
* https://www.nas.nasa.gov/assets/nas/pdf/papers/NAS_Technical_Report_NAS-2015-02.pdf
* https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt
