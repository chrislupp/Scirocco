# Correctness and gradient review — 2026-09-21

Before the fix, the existing tests passed, but gradients were incorrect
across the public interface. The most consequential verified defect is a mismatch between the
geometry used by the viscous primal and the geometry differentiated by its
adjoint. A leading-edge drag derivative is wrong by 18.95% on a converged,
unit-chord NACA 0012 case.

The findings below describe the original implementation before the fixes on
`fix/mfoil-gradient-correctness`. Finding 5 corrects the initial normalization
claim using XFOIL's documented convention. Source references below identify
the reviewed implementation. Temporary reproductions and logs are in
`/tmp/scirocco-review-artifacts/`; the independent Release build is in
`/tmp/scirocco-review-build/`.

## Upstream provenance correction

A subsequent check against [Fidkowski's published Python source](https://websites.umich.edu/~kfid/codes/mfoil/mfoil.py)
(version 2023-06-28) found that the bundled file already contained changes
before this review. Earlier references here to an "unmodified" or "original"
Python primal mean the bundled module without the complex-step monkey patches;
they must not be read as a claim that the bundle was pristine upstream.

The website source was retrieved through the browser because direct downloads
encountered a Cloudflare challenge. The complete source was compared after
normalizing blank lines and trailing whitespace, with an additional Python AST
comparison. This is a source-content comparison, not a byte-for-byte checksum.
The first import, `4ff9e63` (2026-01-13), matches that normalized website source.
All current differences are accounted for by these later commits:

| Change in bundled Python | Commit | Assessment |
|---|---|---|
| `stagpoint_find`: second stagnation-coordinate derivative uses `-G[0]` instead of `-G[1]` | `637233a`, March 21 | Correct derivative of the quotient defining `sstag`; independent of C++. |
| `space_geom`: 10 to 20 iterations, tolerance 1e-6 to 1e-14 | `0aea716`, May 16 | Tighter numerical solve; can change returned wake geometry. |
| `build_glob_sys`: compare `R_x.shape` with the dimensions | `0aea716`, May 16 | Sparse-matrix compatibility fix. |
| `get_cteq`: fix `HK_U` typo to `Hk_U` on the clamp | `0aea716`, May 16 | Correctly zero the derivative of the clamped value. |
| `get_damp`: exponential derivative coefficient -1 to -2 | `0aea716`, May 16 | Correct: derivative of `0.1*exp(-20*h)` is `-2*exp(-20*h)`. |
| `Param.doplot`: True to False | `0aea716`, May 16 | Plotting default only. |
| `residual_transition`: 20 to 30 Newton iterations, tolerance 1e-10 to 1e-14 | `6a56ebb`, May 16 | Changes iteration/stopping behavior, including the iteration-dependent step limit. |
| `Param.rtol`: 1e-10 to 1e-12 | `ff1cb4f`, May 17 | Tighter global Newton convergence criterion. |

The derivative edits follow directly from the Python primal expressions;
the numerical-setting changes must still be distinguished from upstream
behavior. The complete AST comparison found changes only in `Param`,
`stagpoint_find`, `space_geom`, `build_glob_sys`, `residual_transition`,
`get_cteq`, and `get_damp`. In particular, `calc_force`, `spline2d`,
`spline_curvature`, and `TE_info` are identical to the published source.
Thus the force/chord convention and spline construction were not introduced
by historical local edits. The 145 C++ / four Python passing tests reported
below validate the locally corrected/tightened reference, not the untouched
website version. No solver code was changed during this provenance audit.

## Current gradient fixes

The follow-up gradient work on `fix/mfoil-gradient-correctness` addresses the
remaining known derivative defects. The sections after this one retain the
initial review and first-round validation as historical records.

- Panel arc-length derivatives now distinguish the constructor's fixed labels,
  regenerated spline lengths, and cumulative straight-panel lengths. The
  standalone constructor no longer silently differentiates a geometry update
  it did not perform. Philote continues to use the regenerated spline path.
- Inviscid and viscous lift-target sensitivities differentiate the constraint
  `Cl = cltarget`. Output gradients include the resulting angle adjustment;
  target and solved-angle derivatives are returned as well. The inviscid
  public solve now actually enforces its requested lift target.
- Frozen-wake trim retains the wake's tracing angle for geometric derivatives.
  Initial-angle derivatives include wake motion and the subsequent trim
  adjustment. Rebuilt-wake trim instead differentiates the wake at the solved
  angle. This distinction is checked against independently solved Python cases.
- Near-field pressure drag (`cdpi`) now has shape, angle, Mach, Reynolds and
  lift-target derivatives in both solver modes. Total `cd` remains a distinct
  output; its inviscid derivatives remain zero.
- Wake-spacing differentiation uses the same iteration count and tolerance as
  the primal. Newton updates retain Python's stagnation/transition refresh even
  near convergence. Sparse assembly preserves nonzero Jacobian entries rather
  than dropping values below a fixed threshold.
- Cached public gradient accessors reject reads before successful computation
  and after solves, setters, geometry helpers, or mutable state access. A
  retained mutable reference still requires the caller to solve again after
  edits; arbitrary subsequent writes through that reference are not tracked.

The complex-step wrapper was also missing Python's post-update shape-factor
projection, negative shear-state repair, and inviscid/wake rebuild. These were
restored from Python's real `update_state`, not inferred from C++. An independent
unit test compares their real behavior directly with that function. The bundled
`mfoil.py` itself is unchanged by this work.

Seven additional reference cases cover fixed arc labels, negative angle,
cambered viscous flow, inviscid trim, frozen-wake viscous trim, rebuilt-wake
viscous trim, and chord 2. They check all five returned outputs (`Cl`, `Cd`,
`Cm`, `Cdpi`, solved alpha), ten shape directions per case and all four scalar
inputs (initial alpha, Mach, Reynolds, lift target). Combined with the first
reference set this gives 112 shape directions over ten cases. Reynolds
sensitivities are multiplied by Reynolds before comparison.

Unit-chord derivative checks retain `2e-11 + 2e-8*abs(reference)`. The new
nonunit-chord shape checks use `5e-10 + 2e-8*abs(reference)`. Specifically, the
chord-2 node-50 horizontal lift derivative is approximately
`-0.00967396684920` in C++ and `-0.00967396663282` in Python CS, a difference of
`2.16e-10` (relative `2.24e-8`). Tighter nonlinear solves and an independent
40-digit spline calculation did not eliminate it. Its precise origin remains
unresolved; the looser bound is stated explicitly, not a claim of agreement at
the original tighter tolerance. An independent five-point finite-difference
check of the Python primal approached the CS derivative (error `1.4e-8` at
step `1e-5`) before roundoff/solve noise dominated at smaller steps; it cannot
resolve the `2.16e-10` difference. No reference derivative was replaced with a
C++ result.

For chord 2, the reference generator explicitly corrects the upstream Python
pressure integral's repeated in-loop normalization: it integrates with unit
reference chord and then divides lift/pressure drag by chord once and moment
by chord squared. This follows the force-integral definitions. Total drag
retains the existing Python/XFOIL unit-reference convention. This reference
correction does not change the Python or C++ forward implementations. The
nonunit force conventions therefore remain a separate API/theory issue.

These tests establish accuracy for the checked smooth solution branches, not
for arbitrary separated flows, topology/transition changes, singular trim, or
unconverged states. Panel-node sensitivities hold reference chord and moment
reference point fixed and do not differentiate the initial repaneling from
source geometry.

Regenerate the additional references with:

```sh
OPENBLAS_NUM_THREADS=1 python3 scirocco/tests/extract_remaining_gradient_reference.py
```

Final validation of the follow-up: **153/153 C++ tests pass** in the Release
build with Philote, HDF5 and sensitivity-debug tests enabled; **5/5 Python CS
unit tests pass**. Both reference files were regenerated with the corrected
wrapper and record the current Python source hashes. Changed C++ sources pass
clang-format 17 and targeted clang-tidy reports no user-code warnings.

## First-round fixes and validation

`fix/mfoil-gradient-correctness` differentiates the original Python spline
arc-length construction throughout the foil/wake chain. It fixes inviscid
Cd scalar derivatives and rejects invalid primal/adjoint states and stale
interface inputs. Lift-target primal iteration was repaired in the first round; constrained
sensitivities were initially rejected and are now implemented as described above.

The original `mfoil.py` and force-output conventions are unchanged. The CS
wrapper now preserves the trailing-edge absolute-value chain and extends
SciPy's not-a-knot interpolation to complex knot positions. Independent
SciPy comparisons, polynomial reproduction, geometric scaling, and the
unmodified Python primal validate that extension without using C++ answers.

New references contain 42 shape directions across three operating cases,
including LE, stagnation-adjacent and TE coordinates, plus alpha/Mach/Re
checks. C++ end-to-end checks use `2e-11 + 2e-8*abs(reference)` tolerances;
Reynolds derivatives are scaled by Re so absolute tolerance cannot hide them.
Spline arc-length directional derivatives use a `3e-11` absolute tolerance.
Initial curvature-based repaneling is outside the panel-node derivative
contract: call `RefreshPanelGeometry` before shape-sensitive solves, as the
Philote discipline does.

To regenerate and validate the Python references:

```sh
OPENBLAS_NUM_THREADS=1 python3 -m unittest discover -s mfoil -p 'test_*.py'
OPENBLAS_NUM_THREADS=1 python3 scirocco/tests/extract_gradient_contract_reference.py
```

The generator checks convergence, finiteness and agreement with the
unmodified Python real primal before accepting each derivative. The HDF5
file records SHA256 hashes of all three Python source files.

Validation on the fix branch: all 145 C++ tests pass in a fresh Release
build with Philote, HDF5 output, and sensitivity-debug tests enabled. This
includes the interface's actual `Compute`/`ComputePartials` path against the
new Python references. All four independent Python CS tests pass. Changed
C++ files pass clang-format 17; targeted clang-tidy reports no user-code
warnings in the changed core and discipline sources. Python CS validation
is also registered in CI.

## Findings in the original implementation

### 1. [P1] Different arc-length definitions in primal and shape adjoint

Locations: `src/adjoint.cpp:510`, `src/paneling.cpp:124`,
`src/philote/philote_discipline.cpp:145`.

`RefreshPanelGeometry`, called by the optimization discipline, computes
`foil.s` by integrating a cubic spline through the panel nodes. The adjoint
unconditionally uses `DFoilSDFoilXCumulative`, which differentiates a sum of
straight-panel lengths. These are different functions of the coordinates.
The resulting error propagates through the stagnation location, boundary-layer
coordinates, wake spacing, and mass/source coupling.

Reproduction: NACA 0012, 100 panels, alpha=3 degrees, Re=1e6, Ma=0; call
`RefreshPanelGeometry` before solving, as the discipline does. Use mfoil
complex-step with the same spline arc-length construction (ten passes of
five-point quadrature and not-a-knot cubic splines). Node numbers are zero-based.

| Derivative at node 50, vertical coordinate | C++ adjoint | mfoil CS, same spline geometry | Relative error |
|---|---:|---:|---:|
| dCl/dz | -0.002322268644769834 | -0.002255124074985435 | 2.98% |
| dCd/dz | -0.000028948833431194 | -0.000035718116550849 | 18.95% |
| dCm/dz | -0.000528029320219897 | -0.000537651987228819 | 1.79% |

The C++ and complex-step primal coefficients agree within 4e-14 absolute.
Thus the mismatch is not an unconverged or different primal solution.
As a localization experiment, replacing only the relevant `dfoil_s_dx`
columns with complex-step derivatives of the spline calculation gives
dCd/dz=-0.000035718116523625, agreeing within 8e-10 relative. The same
substitution repairs the lift and moment derivatives and the other sampled
nodes. No production source was changed for this experiment.

The current end-to-end shape test explicitly changes the primal to cumulative
arc length (`tests/sensitivity_reference_test.cpp:324`) and skips node 50
(`:359`). It therefore does not test the geometry path used by the discipline.

Required correction: use one arc-length definition in the primal and adjoint.
For fidelity to spline-based mfoil, differentiate the spline construction.
Changing production to cumulative length is an alternative model change that
must be explicit. Validate the actual discipline geometry path, including the
leading edge and trailing-edge-adjacent nodes.

### 2. [P1] Sensitivities are returned after a failed primal solve

Locations: `src/adjoint.cpp:488`, `src/philote/philote_discipline.cpp:156`.

The guard in `ViscousSensitivities` checks that a state matrix exists, but does
not check `glob.conv`. The discipline also discards the boolean returned by
`Solve`. An implicit adjoint evaluated at an unfinished Newton iterate is not
the derivative of the finite-iteration computation, nor a verified derivative
of a converged solution.

Reproduction:

```cpp
Scirocco foil("0012", 100);
foil.GetState().param.niglob = 1;
foil.SetOper(3.0, -1.0, 1e6, 0.0);
bool converged = foil.Solve();              // false; BL residual norm 1.1711
auto sensitivity = foil.ComputeSensitivities(); // succeeds
```

This returns dCl/dalpha=0.1068925716573119 and
dCd/dalpha=0.0007639005159239261 despite the failed solve.

Required correction: reject invalid primal states and propagate a failed
analysis through the discipline instead of exposing ordinary outputs and
partials. The convergence criterion should cover the coupled residual,
including edge-velocity matching and any lift constraint; currently
`SolveCoupled` checks only the three boundary-layer residuals.

### 3. [P1] The discipline substitutes stale or zero gradients on failure

Location: `src/philote/philote_discipline.cpp:178-195`.

If an adjoint throws, `ComputePartials` returns cached sensitivities from a
previous input point, or zeros if there is no cache. Neither is generally the
Jacobian of the current outputs. The warning is conditional on verbosity,
which is zero by default. This can make an optimizer report false stationarity
or take a step based on another design's gradient.

This finding follows directly from the exception path; the Philote server
was not built or run in this review. Required correction: propagate an
explicit derivative/analysis failure. If approximate gradients are wanted,
they need an explicit, observable approximation contract.

Related numerical limitation: `StoreAdjointFactorization` retries a singular
Jacobian with a diagonal ridge. That changes the adjoint equations. The
documentation's claim that gradient bias is simply O(epsilon*norm(R_V)) is
not a valid bound near singularity: the inverse Jacobian's conditioning and
the objective/residual sensitivities also enter. Check the unregularized
adjoint residual and expose when the regularized result is being used.

### 4. [P1] Lift-target mode assembles an uninitialized constraint and can crash

Locations: `src/coupled_solver.cpp:564-577`, `:622-626`, `:1052`.

The lift-constraint row uses `post.cl`, `post.cl_alpha`, and `post.cl_ue`,
but `CalcForce` is only called after `SolveCoupled` in the viscous solve.
For a new solver these fields are zero/empty, leaving the constraint row
without a useful Jacobian. On factorization failure, `dU.setZero()` does
not allocate the still-empty matrix, and `UpdateState` subsequently indexes it.

Reproduction in the fresh Release build:

```cpp
Scirocco foil("0012", 100);
foil.SetOper(3.0, 0.5, 1e6, 0.0);
foil.Solve();
```

This prints `Warning: SparseLU decomposition failed` and exits with SIGSEGV
(exit status 139). Even after the primal failure is repaired, the adjoint
explicitly discards the extra lift-constraint/alpha unknown, so its derivatives
would be fixed-alpha derivatives rather than derivatives of a trimmed solve.

Required correction: update force values/linearizations during trim iteration,
handle a failed linear solve without updating from an invalid `dU`, and either
implement the constrained adjoint or reject sensitivity requests in this mode.

### 5. Nonunit-chord output convention needs clarification (corrected finding)

The original review called the lack of `/chord` in total drag a confirmed bug.
That was too strong. XFOIL deliberately defines `CD = D/q`, using an implied
unit reference chord. Its `CDCALC` uses the same Squire–Young expression as
Python mfoil, without dividing by actual chord. Conventional `Cd = D/(q*c)`
is a different output; they agree at unit chord.

Sources: [Drela and Youngren's XFOIL manual, Units](https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt)
and [XFOIL CDCALC](https://github.com/RobotLocomotion/xfoil/blob/master/src/xfoil.f#L1149-L1185).
The mfoil paper uses conventional chord-normalized coefficients, which is why
comparing the expression with the paper alone suggested a missing factor.

There remains an inconsistency for nonunit geometry: C++ normalizes Cl and Cdf
with actual chord, but its Cd expression follows XFOIL's unit-reference
convention. Python's total drag has that same definition; Python also divides
its pressure-force accumulators by chord inside the panel loop. The observed
factor-of-two Cd change when doubling chord establishes the implemented
scaling, not by itself a violation of the XFOIL convention.

For the gradient fix, Python's forward output definition is authoritative.
Do not insert a C++-only chord normalization. A future change to nonunit-chord
output conventions would need a deliberate, consistent Python/API decision
and separate validation. Derivatives must match the function actually returned.

### 6. [P2] Inviscid Cd scalar derivatives belong to a different output

Locations: `src/adjoint.cpp:363-364`, `:467`; `src/forces.cpp:154`.

`CalcForce` reports Cd=0 in inviscid mode, consistently with the bundled
mfoil. Shape derivatives of Cd are also zero. But `AlphaSensitivities` and
`MachSensitivities` put derivatives of near-field pressure drag, Cdpi, into
the Cd fields passed to the optimizer.

Reproduction: NACA 2412, 80 panels, alpha=5 degrees, Re=0, Ma=0.3:

```text
Cd = 0
Cdpi = -0.005322085194654508
reported dCd/dalpha = -0.001386156255337185
reported dCd/dMa = -0.02524907413674087
```

The derivative of the implemented constant-zero Cd output must be zero.
Required correction: return zero inviscid Cd scalar derivatives and expose
Cdpi scalar derivatives separately if needed.

### 7. [P2] The complex-step reference drops a trailing-edge derivative

Location: `mfoil/mfoil_cs.py:57`.

`tcp = np.real(cross)` removes the imaginary perturbation. Holding the sign
of `abs(cross)` fixed is appropriate away from a sign change; holding its
value fixed is not. The analytic extension is
`np.sign(np.real(cross))*cross`.

For NACA 2412, 80 panels, alpha=5 degrees, the derivative dCl/dx at node 0 is:

```text
C++ analytic:                   -0.7937942200467308
current complex-step reference: -0.85457948
sign-preserving CS reference:   -0.79379422
```

Six trailing-edge/adjacent coordinate checks agree with C++ after preserving
the imaginary part. Thus this particular mismatch is a reference defect,
not evidence that the C++ trailing-edge derivative should be changed.
The checked-in inviscid shape tests use interior nodes and miss this path.

Required correction: repair the complex-step extension, regenerate affected
references, and include trailing-edge coordinates in validation. Every
reference solve should also assert convergence and agreement of its real
primal state before its imaginary part is accepted as a derivative.

## What passed and what it establishes

Built all core tests, including the optional sensitivity isolation tests:

```sh
cmake -S scirocco -B /tmp/scirocco-review-build \
  -DCMAKE_BUILD_TYPE=Release -DSCIROCCO_BUILD_PHILOTE=OFF \
  -DSCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS=ON -DSCIROCCO_BUILD_PROBES=ON
cmake --build /tmp/scirocco-review-build -j 6
ctest --test-dir /tmp/scirocco-review-build --output-on-failure -j 6
```

All 116 tests passed in 36.24 seconds. This covers the supplied primal,
panel influence, closure, transition/residual, force-partial, geometry-chain,
and end-to-end reference cases. The seven end-to-end sensitivity tests were
also run directly to inspect numerical errors:

| Supplied check | Maximum observed relative error |
|---|---:|
| Inviscid alpha, Cl/Cm | 4.4e-12 |
| Inviscid Mach, Cl/Cm | 6.7e-12 |
| Viscous shape, two tested interior nodes, cumulative geometry | 3.0e-11 |
| Viscous alpha, Cl/Cd/Cm | 1.1e-10 |
| Viscous Mach, Cl/Cd/Cm | 3.7e-11 |
| Viscous Re, Cl/Cd/Cm | 8.4e-11 |

Fresh mfoil complex-step alpha and Reynolds checks for NACA 0012 also agree
within 9e-11 relative, with matching real primal outputs. A fresh Mach check
at Ma=0.3 agrees within 4.2e-11 relative for Cl, Cd, and Cm. The spline-based
checks and repaired trailing-edge reference above were generated afresh;
they do not rely on the stored HDF5 answers.

The discrete-adjoint sign and transpose convention is consistent:
`R_U^T lambda = J_U^T`, then `dJ/dp = J_p - lambda^T R_p`.
The pressure integration, Karman–Tsien chain, and wake-state drag
differentiation have strong agreement on the supplied unit-chord cases.
The shape problem is a missing/mismatched upstream geometry chain, not an
incorrect adjoint sign. The current Mach implementation is analytic; the
older status document describing a finite-difference fallback is stale.

These results do not certify arbitrary geometries, separated/transition-switch
states, trim mode, singular Jacobians, or the Philote runtime. An attempted
additional cambered viscous case (NACA 2412, 80 panels, alpha=4 degrees,
Re=1e6, Ma=0.3) did not converge in the Python complex-step solve; its
derivatives were excluded rather than classified as C++ errors. Complex-step
accuracy requires a valid analytic extension, a consistent primal, and a
converged solve. The reference file alone cannot guarantee those conditions.

The immediate priorities are to align the geometry chain, propagate analysis
and derivative failures, and validate the exact interface used by optimization.
