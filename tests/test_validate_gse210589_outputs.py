import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from validate_gse210589_outputs import main  # noqa: E402,F401


class GSE210589SchemaTests(unittest.TestCase):
    def test_frozen_sample_sheet_has_unique_run_and_replicate_keys(self):
        root = Path(__file__).resolve().parents[1]
        samples = pd.read_csv(root / "config/samples/GSE210589_rin34dpa.csv")
        self.assertEqual(len(samples), 6)
        self.assertEqual(samples.sra_run.nunique(), 6)
        self.assertEqual(
            samples.groupby("condition").biological_replicate.nunique().to_dict(),
            {"WT_34DPA": 3, "rin1_34DPA": 3},
        )


if __name__ == "__main__":
    unittest.main()
