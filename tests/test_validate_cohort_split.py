import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from validate_cohort_split import validate_split  # noqa: E402


def row(cohort: str, group: str, role: str) -> dict[str, str]:
    return {"cohort_id": cohort, "independence_group": group, "analysis_role": role, "admission_status": "eligible"}


class CohortSplitTests(unittest.TestCase):
    def test_disjoint_frozen_external_split_passes(self):
        self.assertEqual(
            validate_split([row("development", "lab_a", "development")], [row("external", "lab_b", "frozen_external_validation")]),
            [],
        )

    def test_shared_independence_group_fails(self):
        errors = validate_split(
            [row("development", "shared_lab", "development")],
            [row("external", "shared_lab", "frozen_external_validation")],
        )
        self.assertTrue(any("independence groups" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
