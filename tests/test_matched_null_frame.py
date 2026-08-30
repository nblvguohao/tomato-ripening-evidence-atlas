import sys
import unittest
from pathlib import Path

import pandas as pd


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from enrich_with_matched_nulls import validate_matching_frame  # noqa: E402


class MatchedNullFrameTests(unittest.TestCase):
    def test_observed_frame_is_exact_signature_ids_in_matching_table(self):
        nulls = pd.DataFrame(
            {
                "draw": [1, 1, 2, 2],
                "signature_gene_id": ["sig1", "sig2", "sig1", "sig2"],
                "matched_null_gene_id": ["c1", "c2", "c3", "c4"],
            }
        )
        matched, sizes = validate_matching_frame({"sig1", "sig2", "sig3"}, nulls)
        self.assertEqual(matched, {"sig1", "sig2"})
        self.assertTrue(sizes.eq(2).all())

    def test_rejects_draw_with_missing_control(self):
        nulls = pd.DataFrame(
            {
                "draw": [1, 1, 2],
                "signature_gene_id": ["sig1", "sig2", "sig1"],
                "matched_null_gene_id": ["c1", "c2", "c3"],
            }
        )
        with self.assertRaises(ValueError):
            validate_matching_frame({"sig1", "sig2"}, nulls)


if __name__ == "__main__":
    unittest.main()
