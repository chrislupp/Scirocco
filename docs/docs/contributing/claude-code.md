---
sidebar_position: 1
---

# Claude Code Configuration

Scirocco ships a `CLAUDE.md` at the repository root that Claude Code reads
automatically when you open the project. It covers the build system, test
suite, code architecture, adjoint structure, and the pre-commit style
requirements.

## Adding machine-specific configuration

`CLAUDE.md` is committed and shared. For anything that is specific to your
machine or personal workflow — local tool paths, preferred build directories,
platform-specific commands — create a `CLAUDE.local.md` file at the repository
root alongside `CLAUDE.md`:

```
ProjectTrysil/
├── CLAUDE.md          ← committed, shared
└── CLAUDE.local.md    ← gitignored, yours only
```

`CLAUDE.local.md` is listed in `.gitignore` and will never be committed.
Claude Code imports it automatically via the `@CLAUDE.local.md` directive at
the bottom of `CLAUDE.md`; if the file does not exist the import is silently
skipped, so there is nothing to set up if you have no local overrides.

## What belongs in CLAUDE.local.md

Anything that is true on your machine but not everyone else's:

- Platform-specific tool invocations (e.g. `xcrun` on macOS for `llvm-cov`)
- Non-default build directory paths
- Personal adjoint debugging notes or probe workflows
- Environment variables or shell aliases needed locally

## Example

```markdown
## Coverage (macOS)

​```bash
cmake -S scirocco -B scirocco/build-cov -DCMAKE_BUILD_TYPE=Release \
      -DSCIROCCO_COVERAGE=ON
cmake --build scirocco/build-cov -j
cd scirocco/build-cov && LLVM_PROFILE_FILE="cov-%p.profraw" ctest
xcrun llvm-profdata merge -sparse cov-*.profraw -o cov.profdata
xcrun llvm-cov report scirocco_tests \
    -object scirocco_forward_ref_tests \
    -object scirocco_sensitivity_ref_tests \
    -object scirocco_ping_tests \
    -instr-profile=cov.profdata \
    $(ls ../src/*.cpp)
​```
```
