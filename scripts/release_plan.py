"""Prepare and check a release. Nothing here pushes or publishes.

  release_plan.py prepare minor            # stable release
  release_plan.py prepare minor rc         # prerelease (alpha, beta or rc)
  release_plan.py check --base origin/main # in CI; labels from RELEASE_LABELS

`prepare` makes every versioned change in the checkout: the CMake version,
the CHANGELOG section and, for a stable release, the docs snapshot.
`check` confirms that the checkout holds exactly the release its labels
describe, relative to `--base` (the main revision it will merge into), and
writes the tag to GITHUB_OUTPUT. A tag that already points at HEAD is
reported, not refused, so a rerun after a partial publication goes on.
"""

import argparse
from datetime import date
import json
import os
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from release_docs import prepare_docs  # noqa: E402
from release_version import calculate  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
REPO = "https://github.com/chrislupp/Scirocco"
PROJECT = re.compile(r"(project\(scirocco\s+VERSION\s+)(\d+\.\d+\.\d+)\b")


def cmake_version(text):
    match = PROJECT.search(text)
    if not match:
        raise ValueError("Cannot read the Scirocco project version")
    return match[2]


def set_cmake_version(text, version):
    updated, count = PROJECT.subn(rf"\g<1>{version}", text)
    if count != 1:
        raise ValueError("Expected exactly one Scirocco project version")
    return updated


def changelog_header(result):
    tag = result["tag_version"]
    return tag[1:] if result["is_prerelease"] == "true" else result["new_version"]


def release_notes(text, header):
    """The body of the CHANGELOG section titled `header`, stripped."""
    match = re.search(
        rf"^## \[{re.escape(header)}\][^\n]*\n(.*?)(?=^## \[|^\[[^\]]+\]: |\Z)",
        text, flags=re.M | re.S,
    )
    if not match:
        raise ValueError(f"CHANGELOG.md has no section for {header}")
    return match[1].strip()


def retitle_changelog(text, header, tag, today):
    """Move the [Unreleased] notes under a new version heading."""
    if text.count("## [Unreleased]\n") != 1:
        raise ValueError("CHANGELOG.md needs exactly one ## [Unreleased] heading")
    if not re.search(r"^- ", release_notes(text, "Unreleased"), flags=re.M):
        raise ValueError("The [Unreleased] section of CHANGELOG.md is empty")
    text = text.replace(
        "## [Unreleased]\n", f"## [Unreleased]\n\n## [{header}] - {today}\n", 1
    )
    text = re.sub(r"^\[Unreleased\]: .*\n?", "", text, flags=re.M)
    return (
        text.rstrip("\n")
        + f"\n[Unreleased]: {REPO}/compare/{tag}...HEAD\n"
        + f"[{header}]: {REPO}/releases/tag/{tag}\n"
    )


def git(root, *args):
    return subprocess.check_output(["git", *args], text=True, cwd=root).strip()


def tags_not_at_head(root):
    tags = git(root, "tag", "--list").split()
    at_head = set(git(root, "tag", "--points-at", "HEAD").split())
    return [t for t in tags if t not in at_head], at_head


def prepare(labels, today=None, root=ROOT, docs=prepare_docs):
    cmake = root / "scirocco/CMakeLists.txt"
    changelog = root / "CHANGELOG.md"
    current = cmake_version(cmake.read_text())
    result = calculate(current, labels, git(root, "tag", "--list").split())
    header = changelog_header(result)
    cmake.write_text(set_cmake_version(cmake.read_text(), result["new_version"]))
    changelog.write_text(retitle_changelog(
        changelog.read_text(), header, result["tag_version"],
        (today or date.today()).isoformat(),
    ))
    docs(result["new_version"], result["is_prerelease"] == "true", root / "docs")
    return result


def check(labels, base, root=ROOT):
    """Raise unless the checkout is the release `labels` make from `base`."""
    base_version = cmake_version(
        git(root, "show", f"{base}:scirocco/CMakeLists.txt"))
    tags, at_head = tags_not_at_head(root)
    result = calculate(base_version, labels, tags)
    version, tag = result["new_version"], result["tag_version"]
    problems = []
    head_version = cmake_version((root / "scirocco/CMakeLists.txt").read_text())
    if head_version != version:
        problems.append(
            f"scirocco/CMakeLists.txt says {head_version}; the labels on "
            f"{base_version} give {version}. Run scripts/prepare_release.sh."
        )
    changelog = (root / "CHANGELOG.md").read_text()
    try:
        if not release_notes(changelog, changelog_header(result)):
            problems.append(f"The CHANGELOG section for {tag} is empty")
    except ValueError as error:
        problems.append(str(error))
    if result["is_prerelease"] == "false":
        versions = json.loads((root / "docs/versions.json").read_text())
        config = (root / "docs/docusaurus.config.ts").read_text()
        if not versions or versions[0] != version:
            problems.append(f"docs/versions.json does not start with {version}")
        if f'lastVersion: "{version}"' not in config:
            problems.append(f"docs lastVersion is not {version}")
    if problems:
        raise ValueError("\n".join(problems))
    result["tag_exists"] = str(tag in at_head).lower()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("bump", choices=("major", "minor", "patch"))
    p.add_argument("qualifier", nargs="?", choices=("alpha", "beta", "rc"))
    c = sub.add_parser("check")
    c.add_argument("--base", required=True)
    n = sub.add_parser("notes")
    n.add_argument("header")
    args = parser.parse_args()

    if args.command == "notes":
        print(release_notes((ROOT / "CHANGELOG.md").read_text(), args.header))
        return
    if args.command == "prepare":
        labels = [args.bump, "prerelease" if args.qualifier else "release"]
        labels += [args.qualifier] if args.qualifier else []
        result = prepare(labels)
    else:
        result = check(json.loads(os.environ["RELEASE_LABELS"]), args.base)
    result["changelog_header"] = changelog_header(result)
    output = "".join(f"{key}={value}\n" for key, value in result.items())
    print(output, end="")
    if os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a") as stream:
            stream.write(output)


if __name__ == "__main__":
    try:
        main()
    except ValueError as error:
        print(f"release_plan: {error}", file=sys.stderr)
        sys.exit(1)
