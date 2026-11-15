import React from "react";
import CodeBlock from "@theme/CodeBlock";
import { useDocsVersion } from "@docusaurus/plugin-content-docs/client";

/**
 * The Docker image tag that matches the docs being read: `next` (built from
 * develop) for the unreleased docs, the release version for a release's
 * snapshot.
 */
export function useImageTag(): string {
    const { version } = useDocsVersion();
    return version === "current" ? "next" : version;
}

/** Starts both Philote servers from the published images, in the background. */
export default function DockerRun() {
    const tag = useImageTag();
    return (
        <CodeBlock language="bash">
            {[
                `docker run -d --rm --name scirocco -p 50051:50051 chrislupp/scirocco:${tag}`,
                `docker run -d --rm --name scirocco-kulfan -p 50052:50052 chrislupp/scirocco-kulfan:${tag}`,
            ].join("\n")}
        </CodeBlock>
    );
}
