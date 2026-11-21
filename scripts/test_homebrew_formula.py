import unittest
from homebrew_formula import bump


OLD_SHA = "a" * 64
NEW_SHA = "b" * 64
FORMULA = f"""class Scirocco < Formula
  desc "Airfoil panel code"
  url "https://github.com/chrislupp/Scirocco/archive/refs/tags/v0.3.0.tar.gz"
  sha256 "{OLD_SHA}"
  head "https://github.com/chrislupp/Scirocco.git", branch: "develop"

  resource "extra" do
    url "https://example.com/extra.tar.gz"
    sha256 "{"c" * 64}"
  end
end
"""
NEW_URL = "https://github.com/chrislupp/Scirocco/archive/refs/tags/v0.4.0.tar.gz"


class HomebrewFormulaTest(unittest.TestCase):
    def test_replaces_only_stable_url_and_sha(self):
        result = bump(FORMULA, NEW_URL, NEW_SHA)
        self.assertEqual(
            result,
            FORMULA.replace("v0.3.0", "v0.4.0").replace(OLD_SHA, NEW_SHA),
        )

    def test_rejects_bad_sha(self):
        for sha in ["", "B" * 64, "b" * 63]:
            with self.subTest(sha=sha), self.assertRaises(ValueError):
                bump(FORMULA, NEW_URL, sha)

    def test_rejects_missing_or_repeated_lines(self):
        for formula in [
            FORMULA.replace(f'  sha256 "{OLD_SHA}"\n', ""),
            FORMULA.replace("  url ", "  mirror "),
            FORMULA + f'  sha256 "{OLD_SHA}"\n',
        ]:
            with self.assertRaises(ValueError):
                bump(formula, NEW_URL, NEW_SHA)


if __name__ == "__main__":
    unittest.main()
