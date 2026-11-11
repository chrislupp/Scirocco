import type { Config } from "@docusaurus/types";
import type * as Preset from "@docusaurus/preset-classic";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import { darkCodeTheme, lightCodeTheme } from "./src/prism";

import versions from "./versions.json";

const GITHUB = "https://github.com/chrislupp/Scirocco";

const config: Config = {
    title: "Scirocco",
    tagline: "Viscous airfoil analysis with exact gradients, in C++",
    favicon: "img/favicon.ico",

    future: {
        v4: true,
    },

    url: "https://chrislupp.github.io",
    baseUrl: "/Scirocco/",

    organizationName: "chrislupp",
    projectName: "Scirocco",

    onBrokenLinks: "throw",
    markdown: {
        hooks: {
            onBrokenMarkdownLinks: "throw",
        },
    },

    i18n: {
        defaultLocale: "en",
        locales: ["en"],
    },

    headTags: [
        { tagName: "link", attributes: { rel: "icon", type: "image/svg+xml", href: "/Scirocco/img/favicon.svg" } },
        { tagName: "meta", attributes: { name: "theme-color", content: "#01050f" } },
    ],

    presets: [
        [
            "classic",
            {
                docs: {
                    sidebarPath: "./sidebars.ts",
                    remarkPlugins: [remarkMath],
                    // Emit only KaTeX's visual HTML. The default also emits a
                    // MathML copy, which caused equations to appear twice in
                    // browsers and text extraction when its CSS was blocked.
                    rehypePlugins: [[rehypeKatex, { output: "html" }]],
                    editUrl: `${GITHUB}/tree/develop/docs/`,
                    lastVersion: "current",
                    versions: {
                        current: {
                            label: "Next",
                        },
                    },
                },
                blog: false,
                theme: {
                    // Fonts and KaTeX are bundled locally: no third-party
                    // stylesheet requests, nothing render-blocking.
                    customCss: [
                        "./node_modules/@fontsource-variable/geist/index.css",
                        "./node_modules/@fontsource-variable/geist-mono/index.css",
                        "./node_modules/katex/dist/katex.min.css",
                        "./src/css/custom.css",
                    ],
                },
            } satisfies Preset.Options,
        ],
    ],

    plugins: [
        "./src/plugins/sources.ts",
        [
            "@docusaurus/plugin-client-redirects",
            {
                // Pages that moved in the documentation rework.
                redirects: [
                    { from: "/docs/tutorials", to: "/docs/python" },
                    { from: "/docs/tutorials/inviscid-naca", to: "/docs/python/inviscid-sweep" },
                    { from: "/docs/tutorials/viscous-analysis", to: "/docs/python/viscous-analysis" },
                    { from: "/docs/tutorials/geometry-mods", to: "/docs/python/geometry" },
                    { from: "/docs/tutorials/custom-airfoil", to: "/docs/python/custom-coordinates" },
                    { from: "/docs/tutorials/sensitivities", to: "/docs/python/sensitivities" },
                    { from: "/docs/tutorials/hdf5-output", to: "/docs/python/hdf5-history" },
                    { from: "/docs/getting-started/docker", to: "/docs/python/docker" },
                    { from: "/docs/optimization/server", to: "/docs/python/servers" },
                    { from: "/docs/theory", to: "/docs/concepts/adjoint" },
                    { from: "/docs/concepts/solver-params", to: "/docs/api/solver-params" },
                    { from: "/docs/contributing/claude-code", to: "/docs/contributing" },
                ],
            },
        ],
    ],

    themes: [
        [
            "@easyops-cn/docusaurus-search-local",
            {
                hashed: true,
                indexBlog: false,
                docsRouteBasePath: "/docs",
                searchBarShortcutHint: true,
                searchResultLimits: 8,
                explicitSearchResultPath: true,
            },
        ],
    ],

    themeConfig: {
        image: "img/logo-full-dark.svg",
        colorMode: {
            defaultMode: "dark",
            respectPrefersColorScheme: true,
        },
        navbar: {
            logo: {
                alt: "Scirocco",
                src: "img/logo-nav.svg",
                srcDark: "img/logo-nav-dark.svg",
                width: 136,
                height: 30,
            },
            items: [
                { type: "docSidebar", sidebarId: "docsSidebar", position: "left", label: "Docs" },
                { type: "docSidebar", sidebarId: "pythonSidebar", position: "left", label: "Python" },
                { type: "docSidebar", sidebarId: "apiSidebar", position: "left", label: "Reference" },
                { type: "search", position: "right" },
                // Only once a release has been snapshotted.
                ...(versions.length > 0 ? [{ type: "docsVersionDropdown" as const, position: "right" as const }] : []),
                {
                    href: GITHUB,
                    position: "right",
                    className: "header-github-link",
                    "aria-label": "GitHub repository",
                },
            ],
        },
        footer: {
            links: [
                {
                    title: "Learn",
                    items: [
                        { label: "Quickstart", to: "/docs/getting-started/quickstart" },
                        { label: "Guides", to: "/docs/guides/analysis" },
                        { label: "How it works", to: "/docs/concepts/panel-method" },
                    ],
                },
                {
                    title: "Use",
                    items: [
                        { label: "C++ reference", to: "/docs/api" },
                        { label: "Python and OpenMDAO", to: "/docs/python" },
                        { label: "Troubleshooting", to: "/docs/getting-started/troubleshooting" },
                    ],
                },
                {
                    title: "Project",
                    items: [
                        { label: "GitHub", href: GITHUB },
                        { label: "Changelog", href: `${GITHUB}/blob/main/CHANGELOG.md` },
                        { label: "Contributing", to: "/docs/contributing" },
                    ],
                },
            ],
            copyright: `© ${new Date().getFullYear()} Christopher A. Lupp. Apache License 2.0. A C++ port of mfoil by Krzysztof J. Fidkowski.`,
        },
        prism: {
            theme: lightCodeTheme,
            darkTheme: darkCodeTheme,
            additionalLanguages: ["cmake", "bash", "cpp"],
        },
        tableOfContents: {
            minHeadingLevel: 2,
            maxHeadingLevel: 3,
        },
    } satisfies Preset.ThemeConfig,
};

export default config;
