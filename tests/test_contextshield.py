import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from contextshield import collapse_context_independence_groups, fit_contextshield  # noqa: E402


class ContextShieldTests(unittest.TestCase):
    def _effects(self) -> pd.DataFrame:
        rows = []
        for gene, values in {"core": [.9, .8, .85, .82], "reversed": [.9, .8, -.9, .82]}.items():
            for index, value in enumerate(values):
                rows.append({"cohort_id": f"s{index}", "independence_group": f"g{index}", "gene_id": gene,
                             "signed_effect_rank": value, "n_early": 3, "n_late": 3, "measured": True,
                             "tissue": "fruit", "assay": "RNA-seq", "early_ordinal": 0, "late_ordinal": 3})
        return pd.DataFrame(rows)

    def test_reversal_is_downweighted(self):
        genes, studies = fit_contextshield(collapse_context_independence_groups(self._effects()))
        score = genes.set_index("gene_id")
        self.assertGreater(abs(score.loc["core", "contextshield_score"]), abs(score.loc["reversed", "contextshield_score"]))
        self.assertTrue(np.isfinite(studies["reliability"]).all())
        self.assertAlmostEqual(float(studies["reliability"].sum()), 1.0)

    def test_context_offset_is_reported(self):
        effects = self._effects()
        effects.loc[effects.cohort_id.isin(["s0", "s1"]), "late_ordinal"] = 1
        collapsed = collapse_context_independence_groups(effects)
        genes, _ = fit_contextshield(collapsed)
        self.assertTrue(any(column.startswith("context_offset__") for column in genes.columns))


if __name__ == "__main__":
    unittest.main()
