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
| `release/<tag>` | Release PRs into `main`, made by `scripts/prepare_release.sh` |

1. **Fork or branch** off `develop`.
2. Keep branches focused — one logical change per PR. If a fix requires a
   refactor, split them into separate PRs.
3. **Open a draft PR early** if you want feedback on the approach before
   finishing the implementation.
4. Before marking ready for review:
   - All CI checks pass (build + full test suite).
   - `clang-format` and `clang-tidy` produce no new warnings (see below).
   - `CHANGELOG.md` is updated (see below).
5. Squash-merge is preferred for feature branches; release PRs into `main`
   are merged with a merge commit.

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

`scripts/prepare_release.sh` applies these labels. Release PRs require exactly one `release`/`prerelease` label and one
`major`/`minor`/`patch` label. Prereleases also require exactly one of
`alpha`, `beta`, or `rc`; stable releases must have none.
The bump starts a new version cycle. Once that version has prerelease tags,
subsequent candidates retain it, incrementing the qualifier counter. Moving
from alpha to beta or RC starts that qualifier at 1; moving backwards is rejected.
A `release` PR promotes the active version to stable without another bump.
The bump label applies again after stable publication. Existing stable target
tags are rejected before any files are changed. Test this logic with
`python3 -m unittest discover -s scripts -p 'test_*.py'`.


### Making a release

1. **Prepare.** With `develop` holding everything the release should contain,
   run:

   ```bash
   scripts/prepare_release.sh minor        # or major / patch
   scripts/prepare_release.sh minor rc     # prerelease: alpha / beta / rc
   ```

   It branches `release/<tag>` from `origin/develop`, merges `origin/main`
   in, regenerates the fauXFOIL screenshots (stable releases), and runs
   `scripts/release_plan.py prepare`, which sets the CMake version, moves the
   `[Unreleased]` notes under the new version in `CHANGELOG.md` and, for a
   stable release, snapshots the docs (`docs/versioned_docs/`,
   `docs/versions.json`, `lastVersion`). It builds the docs, checks the
   result, commits, pushes and opens the PR into `main` with its labels.
   It needs `gh`, `npm` and Qt 6 for the screenshots.

2. **Review.** The PR's CI builds and tests the exact commit that ships.
   The Docker workflow builds and smoke-tests both images on amd64 and
   arm64 and publishes them as `candidate-<head commit>`. **Release check**
   (`release.yml`) confirms that the labels, CMake version, CHANGELOG
   section and docs snapshot agree with each other and with `main`, that
   the branch contains `main` (so the merge ships the tested tree), that the
   docs build, and that `RELEASE_TOKEN` can push to this repository and
   `chrislupp/homebrew-tools` and the Docker Hub credentials log in. Merge
   only when every check is green.

3. **Publish.** Merge the PR with a merge commit. `release.yml` rebuilds
   nothing. Its jobs, in order:
   - **Tag**: checks the merge commit has the tested tree and the candidate
     images exist, then pushes a signed tag;
   - **Tag images**: gives the candidate images their version tags (stable
     releases also get the minor, major and `latest` tags);
   - **GitHub release**: creates it from the CHANGELOG section;
   - **Homebrew**: points the `scirocco` and `fauxfoil` formulae at the
     release (stable only);
   - **Merge main into develop**: opens a PR instead if it conflicts;
   - **Check published artifacts** (`release-artifacts.yml`): installs the
     released source archive and exercises the published images through
     the Python Philote client on both architectures.

   Every job skips work that is already done, so a job that fails (a
   network error, say) is fixed by rerunning it. Run
   `release-artifacts.yml` by hand to recheck an existing release.

`scripts/validate_release.sh` runs the whole validation (C++ suites with
Philote, docs, images) on one machine without publishing, for a local check:

```bash
PHILOTE_CPP_DIR=/path/to/Philote-Cpp/build bash scripts/validate_release.sh
```

The ASan and UBSan job runs the default core and Philote suites, including
installed-consumer and RPC checks, in an instrumented build. The Coverage
job runs the default suite in a Clang-instrumented build and uploads the
result to Codecov.

### Main branch protection

`.github/main-protection.json` records the protection of `main`: the
required status checks, an up-to-date branch and resolved conversations,
no force pushes or deletion. It requires no review, because the project has
one maintainer. `scripts/test_main_protection.py` fails if a required check
is not the name of a CI job, so rename a job and its entry together. Apply
the file with:

```bash
gh api --method PUT repos/chrislupp/Scirocco/branches/main/protection \
    --input .github/main-protection.json
```
