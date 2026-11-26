"""Exercise release validation failure paths without builds or publication."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ReleaseValidationTest(unittest.TestCase):
    def run_validation(self, failure="", extra_env=None):
        with tempfile.TemporaryDirectory() as directory:
            temp = Path(directory)
            mock = temp / "command"
            mock.write_text(
                f"#!{sys.executable}\n"
                "import os, sys\n"
                "from pathlib import Path\n"
                "line = Path(sys.argv[0]).name + ' ' + ' '.join(sys.argv[1:])\n"
                "with open(os.environ['COMMAND_LOG'], 'a') as f: f.write(line + '\\n')\n"
                "if os.environ['FAIL_MATCH'] and os.environ['FAIL_MATCH'] in line: sys.exit(17)\n"
                "if line.startswith('docker run '): print('test-container')\n"
                "if line.startswith('docker port '): print('127.0.0.1:12345')\n"
            )
            mock.chmod(0o755)
            for name in ("python3", "cmake", "ctest", "npm", "docker"):
                (temp / name).symlink_to(mock)
            log = temp / "commands"
            result = subprocess.run(
                ["bash", "scripts/validate_release.sh"], cwd=ROOT,
                env={**os.environ, "PATH": f"{temp}:{os.environ['PATH']}",
                     "PHILOTE_CPP_DIR": "/test/philote", "COMMAND_LOG": str(log),
                     "FAIL_MATCH": failure, **(extra_env or {})},
                capture_output=True, text=True,
            )
            return result, log.read_text().splitlines()

    def test_success_validates_both_images_without_publication(self):
        result, commands = self.run_validation()
        self.assertEqual(result.returncode, 0, result.stderr)
        for kind in ("analysis", "geometry"):
            self.assertTrue(any(f"--kind {kind}" in c for c in commands))
        self.assertFalse(any(c.startswith(("docker push", "docker tag")) for c in commands))

    def test_validate_images_zero_skips_docker(self):
        result, commands = self.run_validation(extra_env={"VALIDATE_IMAGES": "0"})
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(any("run build" in c for c in commands))
        self.assertFalse(any(c.startswith("docker ") for c in commands))

    def test_each_failed_gate_stops_validation_and_cleans_up_containers(self):
        for failure in (
            "discover -s scripts", "discover -s mfoil", "cmake -S", "cmake --build",
            "ctest ", "npm --prefix docs ci", "run typecheck", "run build",
            "docker build -f docker/scirocco-server", "--kind analysis",
            "docker build -f docker/kulfan-server", "--kind geometry",
        ):
            with self.subTest(failure=failure):
                result, commands = self.run_validation(failure)
                self.assertNotEqual(result.returncode, 0)
                failed = next(i for i, c in enumerate(commands) if failure in c)
                # Only cleanup may run after a failed validation command.
                self.assertTrue(all(c.startswith(("docker logs", "docker stop"))
                                    for c in commands[failed + 1:]))
                if "--kind" in failure:
                    self.assertIn("docker stop test-container", commands)
