---
title: Contributing
description: How to report bugs and send changes.
---

# Contributing

Bug reports, documentation fixes, test cases, and code are all welcome.

- **Found a bug?** Open an [issue](https://github.com/chrislupp/Scirocco/issues)
  with the commit, compiler, coordinates, operating point, and output with
  `param.verb = 1`.
- **Want to change something bigger?** Open an issue first to talk it through.
- **Ready to send a change?** Open a pull request against `develop`. The full
  guide is
  [CONTRIBUTING.md](https://github.com/chrislupp/Scirocco/blob/develop/CONTRIBUTING.md).

## Before you open a pull request

- CI passes: build and full test suite.
- `clang-format` and `clang-tidy` report nothing new.
- `CHANGELOG.md` has an entry for any feature or fix.
- New solver behaviour comes with a test that fails without it.

Branches are named `feature/<name>`, `fix/<name>`, or `docs/<name>`, with one
logical change per pull request. Feature branches are squash-merged.

## Testing derivatives

Derivative tests compare against complex-step reference data from mfoil,
never against finite differences. Complex step is exact to machine
precision; finite differences aren't. If an adjoint test fails, first check
the primal solve against mfoil, then check each piece of the chain against
its reference before changing any code.

## Working on these docs

The site lives in `docs/` and is built with Docusaurus:

```bash
cd docs
npm ci
npm start
```

Figures and timings on the site come from real runs. After a solver change,
regenerate them with the `docs_figure_data` and `bench` example programs;
`docs/README.md` has the details.

By contributing you agree to license your work under the
[Apache 2.0 License](https://github.com/chrislupp/Scirocco/blob/develop/LICENSE).
