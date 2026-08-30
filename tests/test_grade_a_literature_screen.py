import csv
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCREEN = ROOT / "config/grade_a_literature_screen.csv"


class GradeALiteratureScreenTests(unittest.TestCase):
    def setUp(self):
        with SCREEN.open(newline="", encoding="utf-8") as handle:
            self.rows = list(csv.DictReader(handle))

    def test_all_grade_a_candidates_have_completed_bounded_screen(self):
        self.assertEqual(len(self.rows), 18)
        self.assertEqual(len({row["gene_id"] for row in self.rows}), 18)
        self.assertTrue(all(row["candidate_grade"] == "A" for row in self.rows))
        self.assertTrue(all(not row["literature_screen_status"].startswith("pending") for row in self.rows))

    def test_every_completed_row_has_scope_and_safe_conclusion(self):
        for row in self.rows:
            with self.subTest(gene_id=row["gene_id"]):
                self.assertTrue(row["literature_screen_scope"].strip())
                self.assertTrue(row["literature_screen_conclusion"].strip())
                self.assertEqual(row["screen_date"], "2026-08-27")


if __name__ == "__main__":
    unittest.main()
