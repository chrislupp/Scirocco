"""Calculate a release without mutating the checkout or publishing anything.

A bump label selects a new cycle. While the CMake version has prerelease tags
but no stable tag, its version is retained until stable promotion.
"""

import json
import os
from pathlib import Path
import re
import subprocess


VERSION = r"(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)"
PRERELEASE = re.compile(
    rf"v(?P<version>{VERSION})-(?P<qualifier>alpha|beta|rc)\.(?P<number>[1-9][0-9]*)"
)


def calculate(current, labels, tags):
    if not re.fullmatch(VERSION, current):
        raise ValueError("Invalid current version")
    labels, tags = set(labels), set(tags)
    kinds = labels & {"release", "prerelease"}
    bumps = labels & {"major", "minor", "patch"}
    qualifiers = labels & {"alpha", "beta", "rc"}
    if len(kinds) != 1 or len(bumps) != 1:
        raise ValueError(
            "Exactly one release/prerelease and one major/minor/patch label required"
        )
    prerelease = "prerelease" in kinds
    if len(qualifiers) != int(prerelease):
        raise ValueError(
            "Prerelease requires one qualifier; stable release permits none"
        )

    parsed = [m for tag in tags if (m := PRERELEASE.fullmatch(tag))]
    active = any(m["version"] == current for m in parsed) and f"v{current}" not in tags
    version = current
    if not active:
        parts = list(map(int, current.split(".")))
        index = {"major": 0, "minor": 1, "patch": 2}[next(iter(bumps))]
        parts[index] += 1
        parts[index + 1 :] = [0] * (2 - index)
        version = ".".join(map(str, parts))
    if f"v{version}" in tags:
        raise ValueError(f"Target version v{version} has already been released")

    tag = f"v{version}"
    if prerelease:
        qualifier = next(iter(qualifiers))
        previous = [m for m in parsed if m["version"] == version]
        rank = {"alpha": 0, "beta": 1, "rc": 2}
        if any(rank[m["qualifier"]] > rank[qualifier] for m in previous):
            raise ValueError("Cannot move a prerelease cycle backwards")
        number = 1 + max(
            (int(m["number"]) for m in previous if m["qualifier"] == qualifier),
            default=0,
        )
        tag += f"-{qualifier}.{number}"
    if tag in tags:
        raise ValueError(f"Tag {tag} already exists")
    return {
        "new_version": version,
        "tag_version": tag,
        "is_prerelease": str(prerelease).lower(),
    }


def main():
    cmake = Path("scirocco/CMakeLists.txt").read_text()
    match = re.search(r"project\(scirocco\s+VERSION\s+(\d+\.\d+\.\d+)\b", cmake)
    if not match:
        raise ValueError("Cannot read Scirocco project version")
    labels = json.loads(os.environ["RELEASE_LABELS"])
    tags = subprocess.check_output(["git", "tag", "--list"], text=True).splitlines()
    result = calculate(match[1], labels, tags)
    output = "".join(f"{key}={value}\n" for key, value in result.items())
    print(output, end="")
    with open(os.environ["GITHUB_OUTPUT"], "a") as stream:
        stream.write(output)


if __name__ == "__main__":
    main()
