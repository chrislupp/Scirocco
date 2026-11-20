# Scirocco documentation

The site uses Docusaurus 3, React, TypeScript, and Markdown, with KaTeX for math.
Run these commands from `docs/` with Node.js 20 or newer:

```bash
npm ci
npm start
```

Edit current content in `docs/` (repository path `docs/docs/`), navigation in
`sidebars.ts`, and site settings in `docusaurus.config.ts`. API pages are
handwritten; verify signatures and behavior against the C++ source.

## Writing style

Write for someone who wants to get something done:

- Lead with the result: code, then its output, then the explanation.
- Second person, active voice, short sentences. No "This page explains…"
  preamble and no closing recap.
- State each limitation once, in a "Limits" section on the concept page it
  belongs to, and link to it instead of repeating it.
- One canonical place per piece of setup; other pages link there.
- Every number on the site comes from a real run. Show program output in
  `output` fences, which render as results rather than code.

Custom MDX components, registered in `src/theme/MDXComponents.tsx`:
`<Next>`, `<CpPlot>`, `<PolarPlot>`, `<ConvergencePlot>`, `<PanelPlot>`,
`<OnsetPlot>`, `<ShapeGradientPlot>`, `<TransitionSweepPlot>`, `<E387DragPlot>`,
the verification and validation figures in `src/components/charts/vv.tsx`
(`<CodePolarPlot>`, `<CodeTransitionPlot>`, `<CodeDistributionPlot>`,
`<ExperimentPolarPlot>`, `<GregoryCpPlot>`, `<ColesWadcockCpPlot>`),
`<BenchTable>`, `<QuickstartSource>` (which shows
`scirocco/examples/quickstart.cpp` as built by CI), and `<N2>`.

## Figures and timings

Plots read `src/data/*.json`, written by the `docs_figure_data` example,
except `onset.json`, which comes from
`scirocco/validation/onset_envelope/fit_onset.py --docs`, and
`transition.json`, which comes from
`scirocco/validation/captured_transition/e387_calibrate.py docs`, and
`vv_codes.json` and `vv_experiments.json`, which come from
`scirocco/validation/vv/vv.py` (see its README; it needs XFOIL).
Timings in `src/data/bench.json` come from the `bench` example. Regenerate
both after a change to solver results or performance:

```bash
cmake --build scirocco/build -j --target docs_figure_data bench
scirocco/build/docs_figure_data docs/src/data
scirocco/build/bench        # copy the table into docs/src/data/bench.json
```

Record the machine in `bench.json`. Chart colours are validated for
colour-vision deficiency in both themes; keep series to the `--sc-series-*`
roles in `src/css/custom.css`.

### N2 diagrams

The model diagrams on the Python and optimization pages are
[Prelude](https://chrislupp.github.io/Prelude/) files in `src/data/n2/`. Edit
one in Prelude, or by hand, then render them all into
`src/data/n2/rendered.json` and commit it:

```bash
npm run n2
```

`<N2 name="multipoint" alt="…" />` draws a diagram inline, colored for the
current theme; hovering over or tapping a circle shows the variables it
carries. Keep each diagram in step with its script in `examples/`.

## Checks

Before submitting a change:

```bash
npm run typecheck
npm run build
```

The production build checks internal links and writes `build/`. Draft pages
are excluded from production; publish a page only after checking its content
and adding an appropriate navigation entry. Keep runnable examples and their
documented inputs and expected results together when updating tutorials.

Pushes to `develop` affecting the site trigger `.github/workflows/docs.yml`,
which deploys the built site to GitHub Pages. No manual `gh-pages` branch push
is needed. `versioned_docs/` contains historical snapshots; do not copy current
API claims into an older version without verifying that version's source.

Stable releases automatically snapshot the current docs and sidebar using
Docusaurus `docs:version`, select that version as the default, and commit the
snapshot before release validation. The navbar's version dropdown then lists
it alongside "Next", the develop docs, which move to `/docs/next/`. Prereleases do not create snapshots or
change the default. The release merge back to `develop` triggers deployment;
if that merge conflicts, deployment of the new snapshot waits for resolution.
Tutorials remain runnable examples for users; snapshotting copies documentation
and does not execute tutorials or generate their results.
