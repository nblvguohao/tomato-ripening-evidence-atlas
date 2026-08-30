#!/usr/bin/env python3
"""Create an auditable transcript-to-gene map from Ensembl Plants cDNA FASTA."""
from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    opener = gzip.open if args.fasta.suffix == ".gz" else open
    rows: list[tuple[str, str]] = []
    with opener(args.fasta, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            transcript = line[1:].split()[0]
            match = re.search(r"\bgene:([^\s]+)", line)
            if match is None:
                raise ValueError(f"No gene attribute in FASTA header: {line[:120]}")
            gene = re.sub(r"^gene-", "", match.group(1))
            rows.append((transcript, re.sub(r"\.\d+$", "", gene)))
    if not rows or len({row[0] for row in rows}) != len(rows):
        raise ValueError("Transcript identifiers are empty or non-unique")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("transcript_id\tgene_id\n" + "\n".join("\t".join(row) for row in rows) + "\n", encoding="utf-8")
    print(f"Wrote {len(rows)} transcript-to-gene rows covering {len({row[1] for row in rows})} genes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
