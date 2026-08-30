#!/usr/bin/env python3
"""Download GEO series-matrix archives listed in the manifest.

SRA accessions are deliberately skipped: they require a separate SRA retrieval
and count/metadata processing path, recorded in the protocol before use.
"""

from __future__ import annotations

import argparse
import csv
import urllib.request
from pathlib import Path


def geo_matrix_url(accession: str) -> str:
    prefix = accession[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/matrix/{accession}_series_matrix.txt.gz"


def geo_soft_url(accession: str) -> str:
    prefix = accession[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/soft/{accession}_family.soft.gz"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=60)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with args.manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    for row in rows:
        accession = row["accession"]
        if not accession.startswith("GSE"):
            print(f"SKIP {accession}: requires the SRA/count-matrix workflow.")
            continue
        destination = args.output / f"{accession}_series_matrix.txt.gz"
        soft_destination = args.output / f"{accession}_family.soft.gz"
        if destination.exists() or soft_destination.exists():
            print(f"EXISTS {destination}")
            continue
        url = geo_matrix_url(accession)
        print(f"DOWNLOAD {accession} <- {url}")
        try:
            urllib.request.urlretrieve(url, destination)
        except Exception as exc:  # series matrices are absent for many RNA-seq studies
            destination.unlink(missing_ok=True)
            soft_url = geo_soft_url(accession)
            print(f"NO MATRIX {accession}: {exc}; retrieving SOFT metadata instead.")
            try:
                urllib.request.urlretrieve(soft_url, soft_destination)
            except Exception as soft_exc:
                soft_destination.unlink(missing_ok=True)
                print(f"FAILED {accession}: {soft_exc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
