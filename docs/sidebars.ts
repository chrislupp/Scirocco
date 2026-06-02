import type { SidebarsConfig } from "@docusaurus/plugin-content-docs";

const sidebars: SidebarsConfig = {
    docsSidebar: [
        {
            type: "category",
            label: "Getting Started",
            collapsed: false,
            items: [
                "intro",
                "getting-started/installation",
                "getting-started/quickstart",
            ],
        },
        {
            type: "category",
            label: "Tutorials",
            collapsed: false,
            items: [
                "tutorials/inviscid-naca",
                "tutorials/viscous-analysis",
                "tutorials/custom-airfoil",
                "tutorials/geometry-mods",
                "tutorials/hdf5-output",
                "tutorials/sensitivities",
            ],
        },
        {
            type: "category",
            label: "Concepts",
            collapsed: false,
            items: [
                "concepts/panel-method",
                "concepts/viscous-coupling",
                "concepts/transition",
                "concepts/compressibility",
                "concepts/geometry",
                "concepts/ptc",
                "concepts/solver-params",
            ],
        },
        {
            type: "category",
            label: "Optimization",
            collapsed: false,
            items: [
                "optimization/server",
                "optimization/openmdao-setup",
                "optimization/naca-design",
                "optimization/cst-design",
                "optimization/multipoint",
            ],
        },
        {
            type: "category",
            label: "Theory",
            collapsed: true,
            items: ["theory/index"],
        },
    ],
    apiSidebar: [
        "api/index",
        "api/scirocco",
        "api/foil-state",
        "api/geometry",
        "api/sensitivity",
        "api/hdf5-writer",
        "api/kulfan-discipline",
        "api/internals",
    ],
};

export default sidebars;
