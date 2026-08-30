#!/usr/bin/env python3
"""Prepare GSE78733 with its archived GPL21525 probe annotation.

The series matrix contains numeric Affymetrix transcript-cluster identifiers;
the platform table in the family SOFT supplies the source-bound Solyc mapping.
Only uniquely mapped probe rows are retained before gene-level aggregation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


GENE_PATTERN = re.compile(r"(Solyc\d+g\d+)")


def read_series(path: Path) -> tuple[list[str], pd.DataFrame]:
    titles: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = iter(handle)
        for line in lines:
            if line.startswith("!Sample_title"):
                titles = re.findall(r'"([^"]+)"', line)
            if line.startswith("!series_matrix_table_begin"):
                break
        header = next(lines).rstrip("\n").split("\t")
        rows = []
        for line in lines:
            if line.startswith("!series_matrix_table_end"):
                break
            rows.append(line.rstrip("\n").split("\t"))
    table = pd.DataFrame(rows, columns=[value.strip('"') for value in header])
    table = table.rename(columns={table.columns[0]: "probe_id"})
    for column in table.columns[1:]:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    return titles, table


def platform_mapping(path: Path) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_table = True
                header = next(handle).rstrip("\n").split("\t")
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                values = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                if len(values) == len(header):
                    rows.append(dict(zip(header, values, strict=True)))
    table = pd.DataFrame(rows)
    if not {"ID", "gene_assignment"}.issubset(table.columns):
        raise ValueError("GPL21525 table lacks required ID/gene_assignment columns")
    genes = table["gene_assignment"].fillna("").map(lambda value: sorted(set(GENE_PATTERN.findall(value))))
    table["gene_id"] = genes.map(lambda values: values[0] if len(values) == 1 else pd.NA)
    table["mapping_status"] = genes.map(lambda values: "one_to_one" if len(values) == 1 else ("unmapped" if not values else "ambiguous"))
    return table.rename(columns={"ID": "probe_id"})[["probe_id", "gene_id", "mapping_status"]]


def metadata(sample_ids: list[str], titles: list[str]) -> pd.DataFrame:
    if len(sample_ids) != len(titles):
        raise ValueError("GSE78733 sample titles do not match matrix columns")
    records = []
    for sample_id, title in zip(sample_ids, titles, strict=True):
        lowered = title.lower()
        stage = "breaker" if "_br_" in lowered else ("turning" if "_tu_" in lowered else ("pink" if "_pk" in lowered else "red_ripe"))
        if stage not in {"breaker", "turning", "pink", "red_ripe"}:
            raise ValueError(f"Cannot parse GSE78733 stage from: {title}")
        records.append({"sample_id": sample_id, "study_id": "GSE78733", "title": title, "genotype": "WT", "stage": stage,
                        "stage_ordinal": {"breaker": 1, "turning": 2, "pink": 2, "red_ripe": 3}[stage]})
    frame = pd.DataFrame(records)
    frame["replicate"] = frame.groupby("stage").cumcount() + 1
    return frame


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--series", type=Path, required=True)
    parser.add_argument("--family-soft", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    titles, expression = read_series(args.series)
    sample_ids = expression.columns[1:].tolist()
    mapping = platform_mapping(args.family_soft)
    joined = expression.merge(mapping, on="probe_id", how="left", validate="one_to_one")
    accepted = joined[joined["mapping_status"] == "one_to_one"].dropna(subset=["gene_id"])
    gene_expression = accepted.groupby("gene_id", as_index=False)[sample_ids].median()
    ranks = gene_expression[sample_ids].rank(axis=0, pct=True, method="average")
    output = pd.concat([gene_expression[["gene_id"]], ranks], axis=1)
    sample_metadata = metadata(sample_ids, titles)
    counts = sample_metadata.groupby("stage")["sample_id"].count().to_dict()
    if counts.get("breaker", 0) < 2 or counts.get("red_ripe", 0) < 2:
        raise ValueError("GSE78733 lacks at least two biological samples at the breaker/red-ripe endpoints")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output_dir / "GSE78733_rank01.csv.gz", index=False, compression="gzip")
    sample_metadata.to_csv(args.output_dir / "GSE78733_sample_metadata.csv", index=False)
    audit = {
        "accession": "GSE78733", "platform": "GPL21525", "matrix_probe_rows": int(len(expression)),
        "platform_annotation_rows": int(len(mapping)), "one_to_one_probe_rows": int((mapping.mapping_status == "one_to_one").sum()),
        "accepted_probe_rows": int(len(accepted)), "gene_count": int(len(output)), "sample_count": int(len(sample_ids)),
        "stage_sample_counts": {key: int(value) for key, value in counts.items()}, "mapping_source": str(args.family_soft),
        "admission": "eligible_for_frozen_external_validation",
    }
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    args.audit.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
