import unittest
from release_version import calculate


class ReleaseVersionTest(unittest.TestCase):
    def test_full_cycle_and_next_cycle(self):
        current, tags = "0.1.0", {"v0.1.0"}
        steps = [
            ("alpha", "v0.2.0-alpha.1"),
            ("alpha", "v0.2.0-alpha.2"),
            ("beta", "v0.2.0-beta.1"),
            ("rc", "v0.2.0-rc.1"),
            ("rc", "v0.2.0-rc.2"),
            (None, "v0.2.0"),
            (None, "v0.3.0"),
        ]
        for qualifier, expected in steps:
            labels = (
                ["minor", "prerelease", qualifier]
                if qualifier
                else ["minor", "release"]
            )
            result = calculate(current, labels, tags)
            self.assertEqual(result["tag_version"], expected)
            self.assertEqual(result["is_prerelease"], str(bool(qualifier)).lower())
            current = result["new_version"]
            tags.add(expected)

    def test_bumps(self):
        for bump, expected in [
            ("major", "2.0.0"),
            ("minor", "1.3.0"),
            ("patch", "1.2.4"),
        ]:
            self.assertEqual(
                calculate("1.2.3", ["release", bump], [])["new_version"], expected
            )

    def test_conflicting_or_missing_labels(self):
        for labels in [
            [],
            ["release"],
            ["patch"],
            ["release", "prerelease", "patch"],
            ["release", "major", "minor"],
            ["release", "patch", "rc"],
            ["prerelease", "patch"],
            ["prerelease", "patch", "rc", "beta"],
        ]:
            with self.subTest(labels=labels), self.assertRaises(ValueError):
                calculate("0.1.0", labels, [])

    def test_existing_target_is_rejected(self):
        for kind in ["release", "prerelease"]:
            labels = [kind, "minor"] + (["rc"] if kind == "prerelease" else [])
            with self.assertRaises(ValueError):
                calculate("0.1.0", labels, ["v0.2.0"])

    def test_numeric_candidate_order_and_unrelated_tags(self):
        tags = ["v0.2.0-rc.9", "v0.2.0-rc.10", "v0.3.0-rc.99", "not-a-version"]
        self.assertEqual(
            calculate("0.2.0", ["prerelease", "minor", "rc"], tags)["tag_version"],
            "v0.2.0-rc.11",
        )

    def test_cannot_regress_qualifier(self):
        with self.assertRaises(ValueError):
            calculate("0.2.0", ["prerelease", "minor", "alpha"], ["v0.2.0-rc.1"])


if __name__ == "__main__":
    unittest.main()
