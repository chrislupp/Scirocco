import { themes as prismThemes } from "prism-react-renderer";
import type { Config } from "@docusaurus/types";
import type * as Preset from "@docusaurus/preset-classic";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";

const config: Config = {
    title: "Scirocco",
    tagline: "Inviscid and viscous airfoil analysis in modern C++",
    favicon: "img/favicon.ico",

    future: {
        v4: true,
    },

    url: "https://chrislupp.github.io",
    baseUrl: "/Scirocco/",

    organizationName: "chrislupp",
    projectName: "Scirocco",

    onBrokenLinks: "throw",

    i18n: {
        defaultLocale: "en",
        locales: ["en"],
    },

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
                    editUrl:
                        "https://github.com/chrislupp/Scirocco/tree/develop/docs/",
                    lastVersion: "current",
                    versions: {
                        current: {
                            label: "Next",
                        },
                    },
                },
                blog: false,
                theme: {
                    // Bundle KaTeX locally so visual math styling never
                    // depends on a CDN stylesheet request.
                    customCss: [
                        "./node_modules/katex/dist/katex.min.css",
                        "./src/css/custom.css",
                    ],
                },
            } satisfies Preset.Options,
        ],
    ],

    themeConfig: {
        colorMode: {
            defaultMode: "dark",
            respectPrefersColorScheme: true,
        },
        navbar: {
            title: "Scirocco",
            items: [
                {
                    type: "docSidebar",
                    sidebarId: "docsSidebar",
                    position: "left",
                    label: "Docs",
                },
                {
                    type: "docSidebar",
                    sidebarId: "apiSidebar",
                    position: "left",
                    label: "API",
                },
                {
                    type: "docsVersionDropdown",
                    position: "right",
                },
                {
                    href: "https://github.com/chrislupp/Scirocco",
                    label: "GitHub",
                    position: "right",
                },
            ],
        },
        footer: {
            style: "dark",
            links: [
                {
                    title: "Documentation",
                    items: [
                        {
                            label: "Docs",
                            to: "/docs/intro",
                        },
                        {
                            label: "API Reference",
                            to: "/docs/api",
                        },
                    ],
                },
                {
                    title: "More",
                    items: [
                        {
                            label: "GitHub",
                            href: "https://github.com/chrislupp/Scirocco",
                        },
                        {
                            label: "Changelog",
                            href: "https://github.com/chrislupp/Scirocco/blob/main/CHANGELOG.md",
                        },
                    ],
                },
            ],
            copyright: `Copyright © ${new Date().getFullYear()} Christopher A. Lupp. Built with Docusaurus.`,
        },
        prism: {
            theme: prismThemes.github,
            darkTheme: prismThemes.dracula,
            additionalLanguages: ["cmake", "bash", "cpp"],
        },
    } satisfies Preset.ThemeConfig,
};

export default config;
