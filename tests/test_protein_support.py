import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from process_pxd051570 import eligibility_matched_rate, fisher_greater  # noqa: E402


class ProteinSupportTests(unittest.TestCase):
    def test_fisher_greater_detects_enrichment(self):
        self.assertLess(fisher_greater(20, 5, 10, 100), 1e-8)

    def test_fisher_greater_returns_one_for_no_targets(self):
        self.assertEqual(fisher_greater(0, 10, 0, 100), 1.0)

    def test_matched_rate_uses_only_covariate_eligible_measured_signature(self):
        measured, targets, rate = eligibility_matched_rate(
            {"sig1", "sig2"},
            {"sig1", "sig2", "sig3", "control"},
            {"sig1", "sig3"},
        )
        self.assertEqual((measured, targets), (2, 1))
        self.assertEqual(rate, 0.5)


if __name__ == "__main__":
    unittest.main()
