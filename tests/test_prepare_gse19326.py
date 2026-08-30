import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from prepare_gse19326 import map_signature_probes  # noqa: E402


class GSE19326PreparationTests(unittest.TestCase):
    def test_probe_mapping_requires_six_hits_and_rejects_ambiguous_targets(self):
        common = ["A" * 25, "C" * 25, "G" * 25, "T" * 25, "ACGTA" * 5, "TTGCA" * 5]
        target_one = [
            "AGTCGATCGATCGATCGATCGATCG", "CGTACGTACGTACGTACGTACGTAC",
            "GATCGATCGATCGATCGATCGATCG", "TGCATGCATGCATGCATGCATGCAT",
            "AACCGGTTAACCGGTTAACCGGTTA", "TTGGCCAATTGGCCAATTGGCCAAT",
        ]
        transcripts = {
            "Solyc01g000010": "".join(common + target_one),
            "Solyc01g000020": "".join(common),
        }
        probes = pd.DataFrame({
            "probe_id": ["one"] * 7 + ["ambiguous"] * 6 + ["short"] * 5,
            "sequence": target_one + [target_one[0]] + common + common[:5],
        })
        mapping = map_signature_probes(probes, transcripts).set_index("probe_id")
        self.assertEqual(mapping.loc["one", "mapping_status"], "one_to_one_signature_target")
        self.assertEqual(mapping.loc["one", "gene_id"], "Solyc01g000010")
        self.assertEqual(mapping.loc["ambiguous", "mapping_status"], "ambiguous_within_signature")
        self.assertEqual(mapping.loc["short", "mapping_status"], "unmapped")


if __name__ == "__main__":
    unittest.main()
