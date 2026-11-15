"""Check snapshot selection and failure behavior without publishing a release."""

from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from release_docs import prepare_docs


class ReleaseDocsTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.docs = Path(self.directory.name)
        self.config = self.docs / "docusaurus.config.ts"
        self.original = 'docs: { lastVersion: "current", versions: {} }\n'
        self.config.write_text(self.original)

    @patch("release_docs.subprocess.run")
    def test_stable_release_snapshots_before_selecting_default(self, run):
        def snapshot(*args, **kwargs):
            self.assertEqual(self.config.read_text(), self.original)
        run.side_effect = snapshot
        prepare_docs("0.3.0", False, self.docs)
        run.assert_called_once_with(
            ["npm", "run", "docusaurus", "--", "docs:version", "0.3.0"],
            cwd=self.docs, check=True,
        )
        self.assertEqual(self.config.read_text(),
                         self.original.replace('"current"', '"0.3.0"'))

    @patch("release_docs.subprocess.run")
    def test_prerelease_preserves_existing_default_without_running_npm(self, run):
        self.config.write_text(self.original.replace('"current"', '"0.2.0"'))
        before = self.config.read_text()
        prepare_docs("0.3.0", True, self.docs)
        run.assert_not_called()
        self.assertEqual(self.config.read_text(), before)

    @patch("release_docs.subprocess.run")
    def test_failed_snapshot_does_not_change_default(self, run):
        run.side_effect = subprocess.CalledProcessError(1, "npm")
        with self.assertRaises(subprocess.CalledProcessError):
            prepare_docs("0.3.0", False, self.docs)
        self.assertEqual(self.config.read_text(), self.original)

    @patch("release_docs.subprocess.run")
    def test_invalid_version_or_config_fails_before_snapshot(self, run):
        for version in ("0.3.0-rc.1", "v0.3.0", "../bad"):
            with self.subTest(version=version), self.assertRaises(ValueError):
                prepare_docs(version, False, self.docs)
        for config in ('docs: {}', self.original + self.original):
            self.config.write_text(config)
            with self.subTest(config=config), self.assertRaises(ValueError):
                prepare_docs("0.3.0", False, self.docs)
        run.assert_not_called()
