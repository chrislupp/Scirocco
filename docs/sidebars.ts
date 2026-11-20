import type { SidebarsConfig } from "@docusaurus/plugin-content-docs";

const sidebars: SidebarsConfig = {
    docsSidebar: [
        {
            type: "category",
            label: "Get started",
            collapsible: false,
            items: [
                "intro",
                "getting-started/quickstart",
                "getting-started/installation",
                "getting-started/troubleshooting",
            ],
        },
        {
            type: "category",
            label: "Guides",
            collapsible: false,
            items: ["guides/analysis", "guides/geometry", "guides/gradients", "guides/hdf5"],
        },
        {
            type: "category",
            label: "How it works",
            collapsible: false,
            items: [
                "concepts/panel-method",
                "concepts/viscous-coupling",
                "concepts/transition",
                "concepts/compressibility",
                "concepts/geometry",
                "concepts/ptc",
                "concepts/adjoint",
            ],
        },
        {
            type: "category",
            label: "Verification & validation",
            collapsible: false,
            items: [
                "validation/index",
                "validation/xfoil-mfoil",
                "validation/experiments",
                "validation/transition",
                "validation/performance",
            ],
        },
        "contributing",
    ],
    pythonSidebar: [
        {
            type: "category",
            label: "Python",
            collapsible: false,
            items: ["python/index", "python/docker", "python/servers"],
        },
        {
            type: "category",
            label: "Tutorials",
            collapsible: false,
            items: [
                "python/inviscid-sweep",
                "python/viscous-analysis",
                "python/geometry",
                "python/custom-coordinates",
                "python/sensitivities",
                "python/hdf5-history",
            ],
        },
        {
            type: "category",
            label: "Optimization",
            collapsible: false,
            items: ["optimization/naca-design", "optimization/cst-design", "optimization/multipoint"],
        },
    ],
    apiSidebar: [
        {
            type: "category",
            label: "C++ reference",
            collapsible: false,
            items: [
                "api/index",
                "api/scirocco",
                "api/sensitivity",
                "api/geometry",
                "api/foil-state",
                "api/solver-params",
                "api/hdf5-writer",
                "api/internals",
            ],
        },
        {
            type: "category",
            label: "Conventions",
            collapsible: false,
            items: ["api/conventions", "api/build-options"],
        },
        {
            type: "category",
            label: "Philote",
            collapsible: false,
            items: ["api/kulfan-discipline"],
        },
    ],
};

export default sidebars;
