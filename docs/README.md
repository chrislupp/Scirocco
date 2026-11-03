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
