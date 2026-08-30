import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from studyshield import collapse_independence_groups, random_effects_scores, studyshield_scores  # noqa: E402


class StudyShieldTests(unittest.TestCase):
    def test_shared_study_strata_collapse_to_one_unit(self):
        effects = pd.DataFrame({
            "cohort_id": ["study_A_var1", "study_A_var2", "study_B"],
            "independence_group": ["A", "A", "B"],
            "gene_id": ["g1", "g1", "g1"],
            "signed_effect_rank": [.8, .4, .6],
            "hedges_g": [1.0, .5, .7],
            "hedges_g_variance": [.2, .2, .3],
            "n_early": [3, 3, 3],
            "n_late": [3, 3, 3],
        })
        collapsed = collapse_independence_groups(effects)
        self.assertEqual(collapsed.cohort_id.nunique(), 2)
        shared = collapsed[collapsed.independence_group == "A"].iloc[0]
        self.assertAlmostEqual(shared.signed_effect_rank, .6)
        self.assertEqual(shared.n_early, 6)

    def test_consistent_gene_outranks_heterogeneous_gene(self):
        effects = pd.DataFrame(
            [[.9, .8, .7, .9, .8, .7], [.9, -.8, .7, -.9, .8, -.7]],
            index=["consistent", "heterogeneous"],
            columns=[f"c{i}" for i in range(6)],
        )
        result = studyshield_scores(effects).set_index("gene_id")
        self.assertGreater(abs(result.loc["consistent", "studyshield_score"]), abs(result.loc["heterogeneous", "studyshield_score"]))
        self.assertGreater(result.loc["consistent", "direction_consistency"], result.loc["heterogeneous", "direction_consistency"])

    def test_single_reversal_reduces_score(self):
        base = pd.DataFrame([[.8, .8, .8, .8]], index=["gene"], columns=list("abcd"))
        reversed_one = base.copy()
        reversed_one.iloc[0, -1] = -.8
        score_base = abs(studyshield_scores(base).iloc[0].studyshield_score)
        score_reversed = abs(studyshield_scores(reversed_one).iloc[0].studyshield_score)
        self.assertGreater(score_base, score_reversed)

    def test_missing_cohort_reduces_coverage(self):
        full = pd.DataFrame([[.8, .7, .9]], index=["gene"], columns=list("abc"))
        missing = full.copy()
        missing.iloc[0, -1] = np.nan
        self.assertGreater(studyshield_scores(full).iloc[0].cohort_coverage, studyshield_scores(missing, total_cohorts=3).iloc[0].cohort_coverage)

    def test_random_effects_is_finite(self):
        rows = pd.DataFrame({
            "gene_id": ["g1", "g1", "g1"],
            "hedges_g": [.5, .7, .4],
            "hedges_g_variance": [.2, .2, .3],
        })
        self.assertTrue(np.isfinite(random_effects_scores(rows).loc["g1"]))


if __name__ == "__main__":
    unittest.main()
