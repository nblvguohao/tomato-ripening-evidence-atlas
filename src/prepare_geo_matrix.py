#!/usr/bin/env python3
"""Prepare GEO series matrices in a within-sample rank space for transfer tests.

This deliberately avoids merging absolute microarray values from separate
studies. Only Solyc-labelled rows are retained and gene versions are collapsed.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


def read_series_matrix(path: Path) -> tuple[list[str], pd.DataFrame]:
    titles: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = iter(handle)
        for line in lines:
            if line.startswith("!Sample_title"):
                titles = re.findall(r'"([^"]+)"', line)
            if line.startswith("!series_matrix_table_begin"):
                break
        header = next(lines)
        columns = [value.strip('"\n') for value in header.rstrip("\n").split("\t")]
        rows: list[list[str]] = []
        for line in lines:
            if line.startswith("!series_matrix_table_end"):
                break
            rows.append([value.strip('"\n') for value in line.rstrip("\n").split("\t")])
    frame = pd.DataFrame(rows, columns=columns)
    frame = frame.rename(columns={columns[0]: "gene_id"})
    for column in frame.columns[1:]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return titles, frame


def clean_gene_ids(series: pd.Series) -> pd.Series:
    return series.astype(str).str.extract(r"(Solyc\d+g\d+)", expand=False)


def metadata_for(accession: str, sample_ids: list[str], titles: list[str]) -> pd.DataFrame:
    if len(sample_ids) != len(titles):
        raise ValueError(f"{accession}: sample-title count does not match matrix columns")
    records = []
    for sample_id, title in zip(sample_ids, titles, strict=True):
        lowered = title.lower()
        if accession == "GSE108415":
            genotype = "WT" if "wild type" in lowered else ("LCY1" if "lcy1" in lowered else "LCY3")
            stage = "mature_green" if "mature green" in lowered else ("breaker_plus_10" if "breaker+10" in lowered else "breaker")
            ordinal = {"mature_green": 0, "breaker": 1, "breaker_plus_10": 3}[stage]
            replicate = int(re.search(r"replicate\s+(\d+)", lowered).group(1))
        elif accession == "GSE42783":
            stage = "mature_green" if "_mg_" in lowered else ("turning" if "_t_" in lowered else "red_ripe")
            ordinal = {"mature_green": 0, "turning": 1, "red_ripe": 3}[stage]
            genotype = "WT"
            replicate = int(re.search(r"rep(\d+)", lowered).group(1))
        else:
            raise ValueError(f"Unsupported accession: {accession}")
        records.append({"sample_id": sample_id, "study_id": accession, "title": title, "genotype": genotype,
                        "stage": stage, "stage_ordinal": ordinal, "replicate": replicate})
    return pd.DataFrame(records)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--accession", choices=("GSE108415", "GSE42783"), required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    titles, raw = read_series_matrix(args.input)
    sample_ids = raw.columns[1:].tolist()
    raw["gene_id"] = clean_gene_ids(raw["gene_id"])
    raw = raw.dropna(subset=["gene_id"])
    expression = raw.groupby("gene_id", as_index=False)[sample_ids].median()
    ranks = expression[sample_ids].rank(axis=0, pct=True, method="average")
    ranked = pd.concat([expression[["gene_id"]], ranks], axis=1)
    metadata = metadata_for(args.accession, sample_ids, titles)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(args.output_dir / f"{args.accession}_rank01.csv.gz", index=False, compression="gzip")
    metadata.to_csv(args.output_dir / f"{args.accession}_sample_metadata.csv", index=False)
    print(f"Prepared {args.accession}: {len(ranked):,} Solyc genes and {len(sample_ids)} samples.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
