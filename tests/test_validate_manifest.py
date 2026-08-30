import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from validate_manifest import validate  # noqa: E402


def row(accession: str, group: str, role: str) -> dict[str, str]:
    return {
        "study_id": accession.lower(),
        "accession": accession,
        "source_url": f"https://example.org/{accession}",
        "organism": "Solanum lycopersicum",
        "assay": "microarray",
        "fruit_or_tissue": "fruit",
        "design": "two replicated stages",
        "intervention": "developmental stage",
        "independence_group": group,
        "planned_role": role,
        "status": "eligible",
        "eligibility_note": "frozen before evaluation",
    }


class ManifestRoleTests(unittest.TestCase):
    def test_frozen_external_role_is_allowed(self):
        self.assertEqual(validate([row("GSE1", "external", "frozen_external_test")]), [])

    def test_frozen_external_group_cannot_also_train(self):
        rows = [row("GSE1", "shared", "frozen_external_test"), row("GSE2", "shared", "train")]
        self.assertTrue(any("spans external testing and development" in error for error in validate(rows)))


if __name__ == "__main__":
    unittest.main()
