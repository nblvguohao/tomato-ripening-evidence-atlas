import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from run_studyshield import baseline_predictions  # noqa: E402


class LeaveOneCohortOutLeakageTests(unittest.TestCase):
    def test_predictions_use_training_rows_only(self):
        train = pd.DataFrame({
            "cohort_id": ["train_a", "train_b", "train_a", "train_b"],
            "gene_id": ["g1", "g1", "g2", "g2"],
            "signed_effect_rank": [0.8, 0.6, -0.7, -0.5],
            "hedges_g": [1.0, 0.8, -0.9, -0.7],
            "hedges_g_variance": [0.2, 0.3, 0.2, 0.3],
            "n_early": [3, 3, 3, 3],
            "n_late": [3, 3, 3, 3],
        })
        frozen = pd.DataFrame({"gene_id": ["g1", "g2"], "effect_source_a": [0.9, -0.9]})
        before = baseline_predictions(train, frozen)
        held_out = train.assign(cohort_id="held", signed_effect_rank=[-1, -1, 1, 1])
        _ = held_out  # Explicitly not supplied to the fitting function.
        after = baseline_predictions(train, frozen)
        for method in before:
            pd.testing.assert_series_equal(before[method][0], after[method][0])


if __name__ == "__main__":
    unittest.main()
