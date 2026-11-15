import React from "react";
import CodeBlock from "@theme/CodeBlock";
import { usePluginData } from "@docusaurus/useGlobalData";

/** quickstart.cpp from the C++ tree, without its license header. */
export function useQuickstartCode(): string {
    const { quickstart } = usePluginData("scirocco-sources") as { quickstart: string };
    return quickstart
        .slice(quickstart.indexOf("*/") + 2)
        .replace(/^\s*\/\/ The program shown on the documentation quickstart page\.\n/, "")
        .trim();
}

export default function QuickstartSource({ title = "quickstart.cpp" }: { title?: string }) {
    return (
        <CodeBlock language="cpp" title={title}>
            {useQuickstartCode()}
        </CodeBlock>
    );
}
