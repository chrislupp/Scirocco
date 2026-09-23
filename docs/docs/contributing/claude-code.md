---
sidebar_position: 1
---

# Contributing

Contributions are welcome — bug reports, documentation improvements, new test
cases, and code changes alike.

## Where to start

- **[CONTRIBUTING.md](https://github.com/chrislupp/Scirocco/blob/develop/CONTRIBUTING.md)**
  — the full contributor guide: branching workflow, commit message conventions,
  code style requirements, changelog format, and how to add tests.
- **[Issue tracker](https://github.com/chrislupp/Scirocco/issues)** — report
  bugs, request features, or discuss ideas before writing code.
- **[Pull requests](https://github.com/chrislupp/Scirocco/pulls)** — all PRs
  target the `develop` branch.

## Quick reference

| Branch prefix | Use for |
| --- | --- |
| `feature/<name>` | New features |
| `fix/<name>` | Bug fixes |
| `docs/<name>` | Documentation-only changes |

Before marking a PR ready:
- All CI checks pass (build + full test suite).
- `clang-format` and `clang-tidy` produce no new warnings.
- `CHANGELOG.md` is updated for any feature addition or bug fix.

Squash-merge is preferred for feature branches.

## Adding tests

New solver features must include at least one test that would fail without the
change. Derivative correctness must be validated against mfoil complex-step
(CS) — never finite differences. See the
[CONTRIBUTING guide](https://github.com/chrislupp/Scirocco/blob/develop/CONTRIBUTING.md#adding-tests)
for details.

## License

By contributing you agree that your changes will be licensed under the
[Apache 2.0 License](https://github.com/chrislupp/Scirocco/blob/develop/LICENSE).
