import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from process_gse40257_rin_binding import fisher_two_sided, peak_targets  # noqa: E402


class RINBindingTests(unittest.TestCase):
    def test_dependency_free_fisher_matches_known_table(self):
        odds_ratio, p_value = fisher_two_sided(50, 413, 1347, 33684)
        self.assertAlmostEqual(odds_ratio, 3.0274432826242874)
        self.assertAlmostEqual(p_value, 1.4529161720667275e-10, places=15)

    def test_requires_two_independent_replicates(self):
        import tempfile

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            columns = {"PEAK_FDR": [0.01, 0.10], "FEATURE_ATTR": ["Solyc01g000010", "Solyc01g000020"],
                       "PEAK_SCORE": [2.0, 3.0], "PEAK_ID": [1, 2]}
            first = root / "GSM1_peak.tsv"
            second = root / "GSM2_peak.tsv"
            pd.DataFrame(columns).to_csv(first, sep="\t", index=False)
            pd.DataFrame({**columns, "FEATURE_ATTR": ["Solyc01g000010", "Solyc01g000030"]}).to_csv(second, sep="\t", index=False)
            targets = peak_targets([first, second]).set_index("gene_id")
            self.assertEqual(int(targets.loc["Solyc01g000010", "binding_replicates"]), 2)
            self.assertNotIn("Solyc01g000020", targets.index)


if __name__ == "__main__":
    unittest.main()
