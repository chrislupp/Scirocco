import fs from "fs";
import path from "path";
import type { LoadContext, Plugin } from "@docusaurus/types";

/**
 * Exposes example sources from the C++ tree to the site at build time, so
 * pages show exactly the file that CI compiles.
 */
const FILES = { quickstart: "../scirocco/examples/quickstart.cpp" };

export default function sourcesPlugin(context: LoadContext): Plugin<Record<string, string>> {
    const resolve = (rel: string) => path.resolve(context.siteDir, rel);
    return {
        name: "scirocco-sources",
        getPathsToWatch: () => Object.values(FILES).map(resolve),
        loadContent: async () =>
            Object.fromEntries(
                Object.entries(FILES).map(([key, rel]) => [key, fs.readFileSync(resolve(rel), "utf8")]),
            ),
        contentLoaded: async ({ content, actions }) => actions.setGlobalData(content),
    };
}
