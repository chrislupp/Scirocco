CST multipoint airfoil optimization
==================================

Optimize one airfoil for low drag at Cl = 0.55 and 0.70, with guard conditions
at Cl = 0.30 and 0.80. This example complements the single-point viscous
optimization in ../kulfan_viscous_optimization/.

Running
-------

Use the Python environment with OpenMDAO, modopt/IPOPT, philote_mdo, NumPy,
SciPy, grpc, and Matplotlib installed. requirements.txt shares the optimizer
stack with the existing Kulfan multipoint example. Build the Scirocco Philote
server, then run from the repository root::

    python examples/cst_range_optimization/run.py \
        --binary scirocco/build-philote/src/philote/scirocco_server \
        --output pa_output/cst_multipoint

--binary defaults to the path above; override it for another build directory.
Each run needs a fresh output directory. The example starts its own analysis
servers and closes them when it exits.

The default uses 100 fixed panels, analytic derivatives, and IPOPT with
20 stored BFGS updates, a tolerance of 1e-6, and at most 160 iterations.
It starts from a symmetric NACA 0012 CST fit. A run succeeds only when both
IPOPT and an independent first-order KKT check pass. Failed trial flow solves
are passed to IPOPT for step rejection.

Problem
-------

Re = 1e6, Mach = 0, and Ncrit = 9 are fixed. The objective is::

    minimize  0.5 Cd(Cl = 0.55) + 0.5 Cd(Cl = 0.70)

There are thirteen design variables: nine shape coefficients and four angles
of attack. The operating conditions are:

===========  ================  ==============================
Target Cl    Objective weight  Additional drag constraint
===========  ================  ==============================
0.30         0                 Cd <= 1.05 Cd_reference(0.30)
0.55         0.5               —
0.70         0.5               —
0.80         0                 Cd <= 1.05 Cd_reference(0.80)
===========  ================  ==============================

Lift equals its target and Cm >= -0.08 at each operating condition.

Continuous camber and half-thickness define the normalized geometry::

    h(x) = sqrt(x) (1-x) T(x) + gap*x/2
    m(x) = x (1-x) M(x)
    upper = m+h; lower = m-h

T has five Bernstein coefficients and M has four. The trailing-edge gap is
fixed at 0.00126. The remaining constraints and variable bounds are:

* Leading-edge radius T(0)^2/2 between 0.005 and 0.020 chord.
* Section area >= 80% of the reference area.
* Thickness >= 80% of reference at 21 stations from 0.15 to 0.65 chord,
  and at 0.90 chord.
* Trailing-edge wedge between 4 and 20 degrees.
* Angles between -2 and 10 degrees.
* T coefficients: 0.1 <= T_0 <= 0.2 and 0 <= T_1,...,T_4 <= 0.4.
* M coefficients between -0.4 and 0.4.
* No curvature constraint.

The reference is the continuous NACA 0012 CST fit used by the existing Kulfan
multipoint example. The initial degree-4 thickness distribution is its
least-squares approximation. Reference drag is evaluated at matched lift.

Results and plots
-----------------

Two runs from the symmetric start reproduced the same design, each converging
in 44 iterations and approximately 90 seconds on the development machine:

===========  =============
Target Cl    Optimized Cd
===========  =============
0.30         0.0054063703
0.55         0.0045067560
0.70         0.0047156273
0.80         0.0080185303
===========  =============

The final objective is 0.0046111917. The final geometry and aerodynamic
gradients were checked against central differences.

Plot the initial and optimized geometry and operating-point forces::

    python examples/cst_range_optimization/plot.py \
        pa_output/cst_multipoint/summary.json --output pa_output/cst_multipoint/result

Plot optimality and feasibility in separate figures::

    python examples/cst_range_optimization/plot_convergence.py \
        pa_output/cst_multipoint/summary.json

This writes optimality.png/.pdf, feasibility.png/.pdf, and convergence.csv
beside summary.json. Use --output to select another directory. The figures
use IPOPT's dual infeasibility and logged original-NLP primal infeasibility,
with logarithmic axes.

summary.json contains the initial and final designs, forces, settings,
optimizer history, independent KKT residuals, timings, and source/binary
hashes. Source snapshots, analysis logs, and process manifests are saved in
the output directory.

Other formulations
------------------

For one drag objective at Cl = 0.55 and guards at Cl = 0.30 and 0.80::

    python examples/cst_range_optimization/run.py \
        --formulation single-point-guards --output pa_output/cst_single_point

The default formulation can also be selected explicitly with
--formulation two-point-guards.

For the original eight-point range benchmark::

    python examples/cst_range_optimization/run.py \
        --formulation multipoint --output pa_output/cst_range

Its objective points are Cl = 0.30, 0.40, 0.50, 0.60, 0.70, and 0.80, with
weights 0.08, 0.20, 0.24, 0.2133333333, 0.1866666667, and 0.08. Its guards
are at Cl = 0.25 and 0.85. Geometry constraints and solver settings are shared
by all formulations.

Checks and optional controls
----------------------------

Run the geometry and KKT-checker unit tests without flow servers::

    python examples/cst_range_optimization/test_geometry.py
    python examples/cst_range_optimization/test_optimality.py

The default run checks a mixed shape/angle derivative direction at three
central-difference steps before optimizing. --skip-check skips this diagnostic.
To check every shape coefficient and angle derivative at a saved design::

    python examples/cst_range_optimization/check_gradients.py \
        pa_output/cst_multipoint/summary.json --output pa_output/cst_gradients

Alternative starts are available through --start cambered and --start thicker.
The cambered start adds 0.5% chord peak camber and a smooth +/-3% thickness
coefficient perturbation. The thicker start increases thickness coefficients
by 3%. --seed permits an explicit restart from a saved summary.

Optional function-value validation utilities are in validate.py; the default
run performs one optimization at one fixed panel count.

Background
----------

* https://web.mit.edu/drela/Public/papers/Pros_Cons_Airfoil_Optimization.pdf
* https://www.brendakulfan.com/_files/ugd/169bff_ccc8f39379d444298b9336f8629f49cb.pdf
* https://www.nas.nasa.gov/assets/nas/pdf/papers/NAS_Technical_Report_NAS-2015-02.pdf
