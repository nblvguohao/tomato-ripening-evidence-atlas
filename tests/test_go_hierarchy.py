import sys
import tempfile
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from propagate_go_annotations import parse_obo  # noqa: E402


class GoHierarchyTests(unittest.TestCase):
    def test_parses_is_a_and_part_of(self):
        content = """format-version: 1.2
[Term]
id: GO:0000001
name: parent
namespace: biological_process
[Term]
id: GO:0000002
name: child
namespace: biological_process
is_a: GO:0000001 ! parent
relationship: part_of GO:0000003 ! another
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "go.obo"
            path.write_text(content)
            parents, terms = parse_obo(path)
        self.assertEqual(parents["GO:0000002"], {"GO:0000001", "GO:0000003"})
        self.assertEqual(terms["GO:0000001"][0], "parent")


if __name__ == "__main__":
    unittest.main()
