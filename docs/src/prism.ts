import type { PrismTheme } from "prism-react-renderer";

/** Code colors built from the logo's blues, with one warm hue for literals. */
export const lightCodeTheme: PrismTheme = {
    plain: { color: "#1b2536", backgroundColor: "#f4f7fb" },
    styles: [
        { types: ["comment", "prolog", "doctype", "cdata"], style: { color: "#7a879b", fontStyle: "italic" } },
        { types: ["punctuation"], style: { color: "#6b778c" } },
        { types: ["operator"], style: { color: "#55627a" } },
        { types: ["keyword", "atrule", "important", "selector"], style: { color: "#2c55a8" } },
        { types: ["builtin", "class-name", "namespace", "tag"], style: { color: "#1f6e99" } },
        { types: ["function", "function-variable"], style: { color: "#0b1220", fontWeight: "600" } },
        { types: ["string", "char", "attr-value", "regex", "inserted"], style: { color: "#1d7a5b" } },
        { types: ["number", "boolean", "constant", "symbol"], style: { color: "#a35a12" } },
        { types: ["macro", "directive", "property"], style: { color: "#6a4fb3" } },
        { types: ["variable", "parameter", "attr-name"], style: { color: "#1b2536" } },
        { types: ["deleted"], style: { color: "#b3412a" } },
    ],
};

export const darkCodeTheme: PrismTheme = {
    plain: { color: "#d5deeb", backgroundColor: "#0a1222" },
    styles: [
        { types: ["comment", "prolog", "doctype", "cdata"], style: { color: "#5d6b82", fontStyle: "italic" } },
        { types: ["punctuation"], style: { color: "#7c8aa3" } },
        { types: ["operator"], style: { color: "#9aa6ba" } },
        { types: ["keyword", "atrule", "important", "selector"], style: { color: "#82a9f2" } },
        { types: ["builtin", "class-name", "namespace", "tag"], style: { color: "#b3e0f8" } },
        { types: ["function", "function-variable"], style: { color: "#f4f7fb", fontWeight: "600" } },
        { types: ["string", "char", "attr-value", "regex", "inserted"], style: { color: "#9dd5c0" } },
        { types: ["number", "boolean", "constant", "symbol"], style: { color: "#f0b97a" } },
        { types: ["macro", "directive", "property"], style: { color: "#c3b5f5" } },
        { types: ["variable", "parameter", "attr-name"], style: { color: "#d5deeb" } },
        { types: ["deleted"], style: { color: "#f08a73" } },
    ],
};
