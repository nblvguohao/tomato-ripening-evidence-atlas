import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from build_hr_submission_bundle import annotation_theme, candidate_grade, literature_screen_complete  # noqa: E402


class HRSubmissionBundleTests(unittest.TestCase):
    def test_grade_a_requires_consensus_and_concordant_protein(self):
        row = {
            "all_four_same_direction": "True",
            "protein_support_status": "differential_direction_concordant",
        }
        self.assertEqual(candidate_grade(row)[0], "A")

    def test_partial_replication_is_grade_c_even_with_protein(self):
        row = {
            "all_four_same_direction": "False",
            "protein_support_status": "differential_direction_concordant",
        }
        self.assertEqual(candidate_grade(row)[0], "C")

    def test_known_anchor_theme_precedes_generic_keywords(self):
        row = {
            "gene_id": "Solyc08g005610",
            "ensembl_gene_name": "CYP707A2",
            "gene_description": "Cytochrome P450",
        }
        self.assertEqual(annotation_theme(row)[0], "Hormone and growth signalling")

    def test_screened_bgh2b_is_not_left_with_limited_annotation(self):
        row = {
            "gene_id": "Solyc10g084600",
            "ensembl_gene_name": "",
            "gene_description": "Plant protein 1589 of Uncharacterized protein function",
        }
        self.assertEqual(annotation_theme(row)[0], "Plastid light adaptation")

    def test_generic_annotation_uses_limited_annotation_not_literature_claim(self):
        row = {
            "gene_id": "Solyc00g000000",
            "ensembl_gene_name": "",
            "gene_description": "Uncharacterized protein",
        }
        self.assertEqual(annotation_theme(row)[0], "Limited functional annotation")

    def test_pending_literature_registry_row_is_not_completed_screen(self):
        self.assertFalse(literature_screen_complete("pending_full_alias_and_full_text_screen"))
        self.assertFalse(literature_screen_complete("not_screened"))
        self.assertTrue(literature_screen_complete("direct_functional_study_outside_fruit_ripening"))


if __name__ == "__main__":
    unittest.main()
