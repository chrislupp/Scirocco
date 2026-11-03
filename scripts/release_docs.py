"""Snapshot stable-release documentation before committing and validating it."""

import argparse
from pathlib import Path
import re
import subprocess


def prepare_docs(version, is_prerelease, docs):
    if is_prerelease:
        return
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        raise ValueError("Documentation snapshots require a stable version")

    config = docs / "docusaurus.config.ts"
    original = config.read_text()
    updated, count = re.subn(
        r'lastVersion:\s*"[^"]*"', f'lastVersion: "{version}"', original
    )
    if count != 1:
        raise ValueError("Expected exactly one docs lastVersion setting")

    # Use the installed, lockfile-pinned CLI. Docusaurus refuses duplicate
    # versions instead of replacing an existing historical snapshot.
    subprocess.run(
        ["npm", "run", "docusaurus", "--", "docs:version", version],
        cwd=docs, check=True,
    )
    config.write_text(updated)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("version")
    parser.add_argument("is_prerelease", choices=("true", "false"))
    args = parser.parse_args()
    prepare_docs(
        args.version, args.is_prerelease == "true",
        Path(__file__).resolve().parents[1] / "docs",
    )
