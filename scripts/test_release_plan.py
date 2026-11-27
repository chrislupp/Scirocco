"""Check release preparation and the release check in a scratch repository."""

from datetime import date
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from release_plan import (
    changelog_header, check, cmake_version, prepare, release_notes,
    retitle_changelog, set_cmake_version,
)

CMAKE = "cmake_minimum_required(VERSION 3.18)\nproject(scirocco VERSION {} LANGUAGES CXX)\n"
CHANGELOG = """# Changelog

## [Unreleased]

### Added
- A new thing.

## [0.3.1] - 2026-11-21

### Fixed
- An old thing.

[0.3.1]: https://github.com/chrislupp/Scirocco/releases/tag/v0.3.1
[Unreleased]: https://github.com/chrislupp/Scirocco/compare/v0.3.1...HEAD
"""


def fake_docs(version, is_prerelease, docs):
    """Stand-in for release_docs.prepare_docs without npm."""
    if is_prerelease:
        return
    versions = docs / "versions.json"
    versions.write_text(json.dumps([version] + json.loads(versions.read_text())))
    config = docs / "docusaurus.config.ts"
    config.write_text(f'lastVersion: "{version}"\n')


class TextTest(unittest.TestCase):
    def test_cmake_version_round_trip(self):
        text = CMAKE.format("0.3.1")
        self.assertEqual(cmake_version(text), "0.3.1")
        updated = set_cmake_version(text, "0.4.0")
        self.assertEqual(cmake_version(updated), "0.4.0")
        self.assertIn("cmake_minimum_required(VERSION 3.18)", updated)

    def test_retitle_moves_notes_and_links(self):
        text = retitle_changelog(CHANGELOG, "0.4.0", "v0.4.0", "2026-11-28")
        self.assertEqual(release_notes(text, "0.4.0"), "### Added\n- A new thing.")
        self.assertEqual(release_notes(text, "Unreleased"), "")
        self.assertEqual(release_notes(text, "0.3.1"), "### Fixed\n- An old thing.")
        self.assertEqual(text.count("[Unreleased]: "), 1)
        self.assertTrue(text.endswith(
            "[Unreleased]: https://github.com/chrislupp/Scirocco/compare/v0.4.0...HEAD\n"
            "[0.4.0]: https://github.com/chrislupp/Scirocco/releases/tag/v0.4.0\n"))

    def test_retitle_refuses_empty_unreleased(self):
        empty = retitle_changelog(CHANGELOG, "0.4.0", "v0.4.0", "2026-11-28")
        with self.assertRaisesRegex(ValueError, "empty"):
            retitle_changelog(empty, "0.5.0", "v0.5.0", "2026-11-29")

    def test_prerelease_header_keeps_qualifier(self):
        result = {"new_version": "0.4.0", "tag_version": "v0.4.0-rc.1",
                  "is_prerelease": "true"}
        self.assertEqual(changelog_header(result), "0.4.0-rc.1")


class RepositoryTest(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        (self.root / "scirocco").mkdir()
        (self.root / "docs").mkdir()
        (self.root / "scirocco/CMakeLists.txt").write_text(CMAKE.format("0.3.1"))
        (self.root / "CHANGELOG.md").write_text(CHANGELOG)
        (self.root / "docs/versions.json").write_text('["0.3.1"]')
        (self.root / "docs/docusaurus.config.ts").write_text('lastVersion: "0.3.1"\n')
        self.git("init", "-q", "-b", "main")
        self.git("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qam", "x",
                 "--allow-empty")
        self.git("add", "-A")
        self.git("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "base")
        self.git("tag", "v0.3.1")
        self.git("checkout", "-qb", "release")

    def git(self, *args):
        return subprocess.check_output(["git", *args], cwd=self.root, text=True)

    def commit(self):
        self.git("add", "-A")
        self.git("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "rel")

    def test_prepared_release_passes_check(self):
        result = prepare(["release", "minor"], date(2026, 11, 28), self.root, fake_docs)
        self.assertEqual(result["tag_version"], "v0.4.0")
        self.commit()
        checked = check(["release", "minor"], "main", self.root)
        self.assertEqual(checked["tag_version"], "v0.4.0")
        self.assertEqual(checked["tag_exists"], "false")

    def test_check_refuses_unprepared_branch(self):
        with self.assertRaisesRegex(ValueError, "says 0.3.1.*give 0.4.0"):
            check(["release", "minor"], "main", self.root)

    def test_check_refuses_labels_that_disagree(self):
        prepare(["release", "minor"], date(2026, 11, 28), self.root, fake_docs)
        self.commit()
        with self.assertRaisesRegex(ValueError, "give 0.3.2"):
            check(["release", "patch"], "main", self.root)

    def test_rerun_after_tagging_reports_existing_tag(self):
        prepare(["release", "minor"], date(2026, 11, 28), self.root, fake_docs)
        self.commit()
        self.git("tag", "v0.4.0")
        self.assertEqual(check(["release", "minor"], "main", self.root)["tag_exists"],
                         "true")

    def test_prerelease_skips_docs_snapshot(self):
        result = prepare(["prerelease", "minor", "rc"], date(2026, 11, 28), self.root,
                         fake_docs)
        self.assertEqual(result["tag_version"], "v0.4.0-rc.1")
        self.commit()
        self.assertEqual(json.loads((self.root / "docs/versions.json").read_text()),
                         ["0.3.1"])
        check(["prerelease", "minor", "rc"], "main", self.root)


if __name__ == "__main__":
    unittest.main()
