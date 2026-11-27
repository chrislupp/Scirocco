"""Every status check main requires must be the name of a real CI job.

A renamed job or a changed matrix otherwise leaves main waiting on a check
that never runs, and release PRs cannot merge.
"""

from itertools import product
import json
from pathlib import Path
import re
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[1]


def job_names(job_id, job):
    """Display names of a job, expanding a simple matrix."""
    name = job.get("name", job_id)
    matrix = (job.get("strategy") or {}).get("matrix") or {}
    if not isinstance(matrix, dict):
        return {name}
    keys = [k for k, v in matrix.items() if k != "include" and isinstance(v, list)]
    combos = [dict(zip(keys, values))
              for values in product(*(matrix[k] for k in keys))] if keys else \
        ([] if matrix.get("include") else [{}])
    for extra in matrix.get("include", []):
        matched = [c for c in combos
                   if all(c.get(k, v) == v for k, v in extra.items() if k in keys)]
        if keys and any(k in keys for k in extra):
            for c in matched:
                c.update(extra)
        elif keys:
            for c in combos:
                if all(c.get(k) == v for k, v in extra.items() if k in c):
                    c.update({k: v for k, v in extra.items() if k not in c})
        else:
            combos.append(dict(extra))

    def expand(combo):
        def value(match):
            item = combo.get(match[1])
            return item.get("name") if isinstance(item, dict) else str(item)
        return re.sub(r"\$\{\{\s*matrix\.(\w+)(?:\.name)?\s*\}\}", value, name)

    return {expand(c) for c in combos}


def all_job_names():
    names = set()
    for path in (ROOT / ".github/workflows").glob("*.yml"):
        workflow = yaml.safe_load(path.read_text())
        for job_id, job in (workflow.get("jobs") or {}).items():
            names |= job_names(job_id, job)
    return names


class MainProtectionTest(unittest.TestCase):
    def test_required_checks_are_real_jobs(self):
        config = json.loads((ROOT / ".github/main-protection.json").read_text())
        names = all_job_names()
        missing = [c for c in config["required_status_checks"]["contexts"]
                   if c not in names]
        self.assertEqual(missing, [], f"not a CI job name; jobs are {sorted(names)}")

    def test_no_review_required(self):
        config = json.loads((ROOT / ".github/main-protection.json").read_text())
        self.assertIsNone(config["required_pull_request_reviews"])


if __name__ == "__main__":
    unittest.main()
