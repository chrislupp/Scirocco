# Contributing to Scirocco

Thank you for your interest in contributing. This document covers everything
you need to get a change from idea to merged PR.

## Table of Contents

- [Branching and PR workflow](#branching-and-pr-workflow)
- [Code style](#code-style)
- [Changelog requirement](#changelog-requirement)
- [Testing](#testing)
- [Adjoint regression debugging](#adjoint-regression-debugging)
- [Commit messages](#commit-messages)

---

## Branching and PR workflow

| Branch | Purpose |
|--------|---------|
| `main` | Stable, tagged releases only |
| `develop` | Integration branch — all PRs target this |
| `feature/<name>` | New features |
| `fix/<name>` | Bug fixes |
| `docs/<name>` | Documentation-only changes |

1. **Fork or branch** off `develop`.
2. Keep branches focused — one logical change per PR. If a fix requires a
   refactor, split them into separate PRs.
3. **Open a draft PR early** if you want feedback on the approach before
   finishing the implementation.
4. Before marking ready for review:
   - All CI checks pass (build + full test suite).
   - `clang-format` and `clang-tidy` produce no new warnings (see below).
   - `CHANGELOG.md` is updated (see below).
5. Squash-merge is preferred for feature branches; merge commits are used for
   release integrations from `develop` → `main`.

---

## Code style

The project enforces consistent formatting and static analysis via
`.clang-format` and `.clang-tidy` at the repo root. Run both before every
commit.

### Formatting

```bash
clang-format -i scirocco/src/*.cpp scirocco/include/scirocco/**/*.h
```

No manual formatting decisions — let clang-format own it entirely.

### Static analysis

A compile database must exist first (any configured build directory works):

```bash
cmake -S scirocco -B scirocco/build -DCMAKE_BUILD_TYPE=Release
run-clang-tidy -p scirocco/build \
    -header-filter='scirocco/include/' \
    'scirocco/src/.*\.cpp$'
```

Fix all new warnings before opening a PR. Do not suppress warnings with
`NOLINT` unless there is a documented reason in the same line comment.

### General conventions

- **`#pragma once`** for all header guards (not `#ifndef`).
- **Doxygen-style `///` comments** for public API declarations.
- **C++17** — no C++20 features.
- Avoid raw owning pointers; prefer value types, `std::unique_ptr`, or
  Eigen containers.
- Do not add comments that restate what the code does. Reserve comments for
  non-obvious invariants, hidden constraints, or workarounds for specific
  external bugs.

---

## Changelog requirement

Every PR that adds a feature, changes existing behaviour, or fixes a bug
**must** include a `CHANGELOG.md` entry under `## [Unreleased]`. Pure
documentation or formatting PRs are exempt.

### Format

```markdown
## [Unreleased]

### Added
- `ClassName::method_name` — brief description of what was added and why it
  matters to users. Wrap at 80 columns.

### Changed
- Modified behaviour, renamed symbols, tightened tolerances, default value
  changes.

### Fixed
- Incorrect results, crashes, or stale documentation that misrepresented
  behaviour.
```

One bullet per logical change. Lead with the affected symbol or subsystem in
backticks when it helps readers locate the change.

---

## Testing

### Running the default suite

```bash
cmake -S scirocco -B scirocco/build -DCMAKE_BUILD_TYPE=Release
cmake --build scirocco/build -j
ctest --test-dir scirocco/build --output-on-failure
```

Test count and runtime depend on the enabled features and installed
dependencies. Inspect the configured suite with `ctest --test-dir
scirocco/build -N`. All configured tests must pass before a PR is ready for
review.

With GTest and HDF5 installed, the default suite builds these executables,
plus an installed-consumer check. Without HDF5 the reference tests are
omitted; without GTest only the installed-consumer check is registered.
Configure with `-DSCIROCCO_BUILD_PHILOTE=ON` to add the Philote discipline
and RPC tests.

| Executable | Purpose |
|---|---|
| `scirocco_solver_strategy_tests` | Solver strategies, input validation, and solver reuse |
| `scirocco_tests` | Spline, geometry, and the mfoil primal reference suite |
| `scirocco_forward_ref_tests` | Forward `R_x` / `R_U` partials vs mfoil |
| `scirocco_ping_tests` | Analytic Jacobians vs mfoil-CS references (BL, residuals, BuildGlobSys, WakeSys, forces) |
| `scirocco_sensitivity_ref_tests` | End-to-end adjoint `dJ/d{x, α, Ma, Re}` vs complex-step references |
| `scirocco_gradient_contract_tests` | Gradient API contracts against references |
| `scirocco_remaining_gradient_tests` | Additional and high-precision gradient references |
| `scirocco_viscous_chain_tests` | Geometric derivative chains against complex-step references |

The opt-in sensitivity-debug suites are described under
[Adjoint regression debugging](#adjoint-regression-debugging).

### Running a single test case

```bash
# One executable
scirocco/build/scirocco_ping_tests

# One test case within an executable
scirocco/build/scirocco_sensitivity_ref_tests \
    --gtest_filter="SensitivityRefTest.AdjointAlpha"
```

### Memory and undefined-behavior checks

```bash
cmake -S scirocco -B scirocco/build-sanitized -DCMAKE_BUILD_TYPE=Debug \
    -DSCIROCCO_SANITIZERS=ON -DSCIROCCO_NATIVE_ARCH=OFF
cmake --build scirocco/build-sanitized -j2
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
    ctest --test-dir scirocco/build-sanitized --output-on-failure
```

The sanitizer CI job also enables Philote and HDF5 output to exercise RPC and
writer paths. Sanitized builds are diagnostic builds, not release artifacts.

### Coverage measurement

Coverage uses Clang compiler flags directly; there is no `SCIROCCO_COVERAGE`
CMake option. CI runs the same measurement in the `Coverage` job and uploads
it to Codecov. To reproduce it locally with Clang and matching LLVM tools on
PATH, and GTest and HDF5 installed:

```bash
CC=clang CXX=clang++ cmake -S scirocco -B scirocco/build-cov \
    -DCMAKE_BUILD_TYPE=Release -DSCIROCCO_NATIVE_ARCH=OFF \
    -DCMAKE_CXX_FLAGS="-fprofile-instr-generate -fcoverage-mapping \
        -fprofile-list=$PWD/.github/coverage-profile-list.txt"
cmake --build scirocco/build-cov -j2
(
    cd scirocco/build-cov
    LLVM_PROFILE_FILE="cov-%p.profraw" ctest --output-on-failure -E '^InstalledConsumer$'
    llvm-profdata merge -sparse cov-*.profraw -o cov.profdata
    llvm-cov report scirocco_tests \
        -object scirocco_solver_strategy_tests \
        -object scirocco_forward_ref_tests \
        -object scirocco_sensitivity_ref_tests \
        -object scirocco_gradient_contract_tests \
        -object scirocco_remaining_gradient_tests \
        -object scirocco_viscous_chain_tests \
        -object scirocco_ping_tests \
        -instr-profile=cov.profdata ../src/*.cpp
)
```

The profile list leaves Eigen, the standard library, GoogleTest, gRPC, and
protobuf uninstrumented. Their inlined counters made the tests about 17x
slower, and the report excludes those headers anyway.

The installed-consumer check is excluded because its separate build does not
inherit the global coverage flags. Use `xcrun llvm-profdata` and
`xcrun llvm-cov` with Apple's bundled Clang tools.

### Adding tests

- New solver features must be accompanied by at least one test that would fail
  without the change.
- Derivative correctness must be validated against **mfoil complex-step (CS)**
  references — never finite differences. CS is exact to machine precision; FD
  is not an acceptable reference.
- If you add a new HDF5 reference dataset, document how to regenerate it in a
  comment near the test.

---

## Adjoint regression debugging

When `scirocco_sensitivity_ref_tests` fails, use the layered debug build
before modifying any adjoint code.

### Step 1 — build the debug suite

```bash
cmake -S scirocco -B scirocco/build-debug -DCMAKE_BUILD_TYPE=Debug \
      -DSCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS=ON \
      -DSCIROCCO_BUILD_PROBES=ON
cmake --build scirocco/build-debug -j
```

### Step 2 — isolate the broken phase

Run each chain-piece executable in order:

| Executable | Adjoint phase |
|------------|--------------|
| `scirocco_viscous_chain_tests` | Forward geometric chain: d(foil.s)/d(foil.x) |
| `scirocco_sensitivity_tests` | Phase 1 — influence derivatives: d(AIC, cp, circ)/d(foil.x) |
| `scirocco_viscous_residual_tests` | Phase 2 — residual partials: d(residual)/d(state) |
| `scirocco_viscous_output_ping_tests` | Phase 3 — output partials: d(Cl, Cd, Cm)/d(state) |

The first executable that fails is the phase to investigate.

### Step 3 — use probe binaries for ad-hoc inspection

The 14 `probe_*` binaries (not part of `ctest`) are standalone diagnostic
programs for comparing individual adjoint terms against CS. Run them directly
from `scirocco/build-debug/`.

### Step 4 — validate the primal first

Before assuming an adjoint error, verify the primal solve matches mfoil. A
primal discrepancy will cascade through every adjoint phase.

---

## Commit messages

Follow the existing style visible in `git log`:

```
Category: short imperative summary (≤72 chars)

Optional body explaining why, not what. Wrap at 72 columns.
```

Common categories: `Fix`, `Add`, `Refactor`, `Docs`, `CI`, `Test`.

Do not reference issue numbers in the subject line; put them in the body if
needed.

## Release version progression

Release PRs require exactly one `release`/`prerelease` label and one
`major`/`minor`/`patch` label. Prereleases also require exactly one of
`alpha`, `beta`, or `rc`; stable releases must have none.
The bump starts a new version cycle. Once that version has prerelease tags,
subsequent candidates retain it, incrementing the qualifier counter. Moving
from alpha to beta or RC starts that qualifier at 1; moving backwards is rejected.
A `release` PR promotes the active version to stable without another bump.
The bump label applies again after stable publication. Existing stable target
tags are rejected before any files are changed. Test this logic with
`python3 -m unittest discover -s scripts -p 'test_*.py'`.


### Release validation and publication

The release workflow pins the merged PR commit and requires it to remain the tip
of `main`. Stable releases also run Docusaurus `docs:version`, set the snapshot
as the default docs version, and include the generated docs, sidebar, version
list, and config in the release commit. Prereleases leave snapshots and the
default docs version unchanged. An existing snapshot is never overwritten;
snapshot-generation failures stop the release before publication.

It prepares and commits these changes locally, then runs
`scripts/validate_release.sh` against that prepared revision. This gate runs the
Python tests, full C++ sensitivity and RPC suite, installed-consumer check,
and documentation build/typecheck. A second job builds both Docker images on
native `linux/amd64` and `linux/arm64` runners from the prepared commit, runs
the malformed-request recovery and valid-request smoke tests against each, and
pushes each tested image to Docker Hub by digest only, without a tag. Any
failure stops publication.

Only after every image passes does the workflow push the release commit/tag,
tag the multi-arch manifests built from those digests, and create the GitHub
release. The atomic Git
push uses an explicit lease on the original `main` revision to reject concurrent
changes. Stable releases receive version, minor, major, and `latest` image
tags; prereleases receive only their full version tag. Manual tag pushes do not
publish images. The Docker workflow builds and smoke-tests both images on
amd64 and arm64 runners on PRs and branch pushes. On pushes to `develop` it
pushes each tested image by digest and, once all four have passed, tags the
two multi-arch manifests as `next`. No other branch publishes images.

For a local validation run without publication, install the dependencies listed
in `scripts/validate_release.sh`, including a running Docker daemon, then run:

```bash
PHILOTE_CPP_DIR=/path/to/Philote-Cpp/build bash scripts/validate_release.sh
```

The workflow uses the existing `RELEASE_TOKEN`, `DOCKERHUB_USERNAME`, and
`DOCKERHUB_TOKEN` secrets. Repository protection must permit the release token to
push the validated version commit. GitHub and Docker Hub publication cannot be
one transaction: if publication fails after validation, inspect the existing tag
and image tags before retrying. Never rerun version calculation blindly after a
partially published release.

## Release qualification

Before promoting `develop` to a stable release:

1. Merge the intended changes into `develop` after their checks pass.
2. Open a `develop` → `main` prerelease PR with one bump label (`major`, `minor`,
   or `patch`), `prerelease`, and `rc`. Have it reviewed and require all checks.
3. The release workflow validates the prepared commit and images before
   publishing. The **Published release artifacts** workflow then downloads the
   tagged source archive, tests a relocated CMake install, pulls both full-version
   Docker images, and exercises them through the Python Philote client.
4. Require that artifact workflow to pass before opening a stable-release PR.
   Its `workflow_dispatch` input allows an existing release tag to be rechecked;
   it does not publish anything. Inspect recorded image digests when diagnosing
   a failed run. A post-publication failure must block promotion, even though it
   cannot undo an already-published prerelease.

The **ASan and UBSan** job runs the default core and Philote suites, including
installed-consumer and RPC checks, in an instrumented build. Release builds
remain uninstrumented. The **Coverage** job runs the default suite in a
Clang-instrumented build and uploads the result to Codecov.

### Main branch protection

`.github/main-protection.json` records the protection configuration. It requires
one approving review, dismissal of stale reviews, resolved conversations, an
up-to-date branch, and the listed CI checks. Force pushes and deletion are
blocked for ordinary contributors. An administrator can apply it with:

```bash
gh api --method PUT repos/chrislupp/Scirocco/branches/main/protection \
    --input .github/main-protection.json
```

Administrator enforcement is deliberately disabled because `RELEASE_TOKEN` is
owned by the repository owner, `chrislupp`, and the release workflow pushes the
validated version commit directly. This exception also permits administrators'
manual pushes; a personal token cannot be distinguished from its owner's other
Git operations. Personal repositories do not support per-user PR bypass lists.
For an automation-only exception, migrate release authentication to a dedicated
GitHub App and use rulesets with that App as the bypass actor. Never broaden the
exception to ordinary contributors or disable checks to publish a release.
