import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from summarize_gse210589_evidence import enrichment_table  # noqa: E402


class GSE210589EvidenceSummaryTests(unittest.TestCase):
    def test_only_admitted_and_tested_targets_enter_enrichment(self):
        differential = pd.DataFrame({
            "gene_id": ["sig_de", "sig_no", "other_de", "other_no"],
            "differential": [True, False, True, False],
        })
        targets = pd.DataFrame({
            "gene_id": ["sig_de", "sig_no", "other_de", "other_no", "rejected"],
            "admitted_binding": [True, True, True, True, False],
        })
        signature = pd.DataFrame({"gene_id": ["sig_de", "sig_no", "rejected"]})
        table, tested = enrichment_table(differential, targets, signature)
        self.assertEqual(table, [[1, 1], [1, 1]])
        self.assertEqual(tested, 4)


if __name__ == "__main__":
    unittest.main()
