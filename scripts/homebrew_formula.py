"""Point the Homebrew formula at a release's source archive.

Usage: python3 scripts/homebrew_formula.py <tag> <path/to/scirocco.rb>

The release workflow runs this against a checkout of chrislupp/homebrew-tools
after it publishes a stable release. It downloads the tag's GitHub source
archive, hashes it, and rewrites the formula's stable `url` and `sha256`.
"""

import hashlib
from pathlib import Path
import re
import sys
import time
import urllib.request


ARCHIVE = "https://github.com/chrislupp/Scirocco/archive/refs/tags/{tag}.tar.gz"
URL_LINE = re.compile(r'^  url "[^"]*"$', re.MULTILINE)
SHA_LINE = re.compile(r'^  sha256 "[0-9a-f]{64}"$', re.MULTILINE)


def bump(formula, url, sha256):
    """Return the formula text with its stable url and sha256 replaced."""
    if not re.fullmatch(r"[0-9a-f]{64}", sha256):
        raise ValueError(f"Invalid sha256: {sha256!r}")
    for name, pattern in [("url", URL_LINE), ("sha256", SHA_LINE)]:
        if len(pattern.findall(formula)) != 1:
            raise ValueError(f"Expected exactly one top-level {name} line")
    formula = URL_LINE.sub(f'  url "{url}"', formula)
    return SHA_LINE.sub(f'  sha256 "{sha256}"', formula)


def archive_sha256(url, attempts=5):
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                return hashlib.sha256(response.read()).hexdigest()
        except OSError:
            if attempt == attempts:
                raise
            time.sleep(10 * attempt)


def main(tag, path):
    url = ARCHIVE.format(tag=tag)
    sha256 = archive_sha256(url)
    formula = Path(path)
    formula.write_text(bump(formula.read_text(), url, sha256))
    print(f"{formula}: {url} sha256 {sha256}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__.strip().splitlines()[2])
    main(sys.argv[1], sys.argv[2])
