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

    stylesheets: [
        {
            href: "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css",
            type: "text/css",
            integrity:
                "sha384-nB0miv6/jRmo5RLHO8BIp/8hwC1slNFDuv3VUgI4A7O/lEKmpRgiuNOQI2bRpFB",
            crossorigin: "anonymous",
        },
    ],

    presets: [
        [
            "classic",
            {
                docs: {
                    sidebarPath: "./sidebars.ts",
                    remarkPlugins: [remarkMath],
                    rehypePlugins: [rehypeKatex],
                    editUrl:
                        "https://github.com/chrislupp/Scirocco/tree/develop/docs/",
                    lastVersion: "0.2.0",
                    versions: {
                        current: {
                            label: "Next",
                        },
                    },
                },
                blog: false,
                theme: {
                    customCss: "./src/css/custom.css",
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
