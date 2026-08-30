import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from matched_crosslayer_enrichment import exact_conditional_odds_ratio_ci, raw_table  # noqa: E402


class ExactConditionalIntervalTests(unittest.TestCase):
    def test_interval_contains_cross_product_odds_ratio(self):
        lower, upper = exact_conditional_odds_ratio_ci(50, 413, 1347, 33684)
        self.assertLess(lower, 3.0274432826242874)
        self.assertGreater(upper, 3.0274432826242874)

    def test_raw_table_records_exact_interval_and_frozen_counts(self):
        table = raw_table({"a", "b", "c"}, {"a", "b", "c", "d", "e", "f"}, {"a", "d", "e"})
        self.assertEqual([table[key] for key in ("signature_target", "signature_not_target", "nonsignature_target", "nonsignature_not_target")], [1, 2, 2, 1])
        self.assertLess(table["exact_conditional_95_ci"][0], table["odds_ratio"])
        self.assertGreater(table["exact_conditional_95_ci"][1], table["odds_ratio"])


if __name__ == "__main__":
    unittest.main()
