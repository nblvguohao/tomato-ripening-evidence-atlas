import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from prepare_gse78733 import platform_mapping  # noqa: E402


class GSE78733PreparationTests(unittest.TestCase):
    def test_platform_mapping_keeps_only_unique_solyc_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "family.soft.gz"
            content = """^PLATFORM = GPL21525
!platform_table_begin
ID\tgene_assignment
probe1\tSolyc01g000010.1.1
probe2\tSolyc01g000020.1.1 /// Solyc01g000030.1.1
probe3\t---
!platform_table_end
"""
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(content)
            mapped = platform_mapping(path).set_index("probe_id")
            self.assertEqual(mapped.loc["probe1", "gene_id"], "Solyc01g000010")
            self.assertEqual(mapped.loc["probe2", "mapping_status"], "ambiguous")
            self.assertEqual(mapped.loc["probe3", "mapping_status"], "unmapped")


if __name__ == "__main__":
    unittest.main()
