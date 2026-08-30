#!/usr/bin/env python3
"""Prepare the GSE128739 pericarp count table without treating lanes as replicates."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


STAGE_ORDINAL = {8: -5, 15: -4, 21: -3, 28: -2, 34: -1, 41: 0, 49: 1, 50: 2, 53: 3}


def core_gene(value: object) -> str | None:
    match = re.search(r"(Solyc\d+g\d+)", str(value))
    return match.group(1) if match else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    raw = pd.read_excel(args.input, header=None)
    stage_values = pd.to_numeric(raw.iloc[0, 1:], errors="raise").astype(int).tolist()
    biological_ids = pd.to_numeric(raw.iloc[1, 1:], errors="raise").astype(int).tolist()
    lane_ids = raw.iloc[2, 1:].astype(str).tolist()
    if len(stage_values) != len(biological_ids) or len(lane_ids) != len(stage_values):
        raise ValueError("Header rows have inconsistent lane counts.")

    counts = raw.iloc[3:, :].copy()
    counts.columns = ["source_gene", *range(len(stage_values))]
    counts["gene_id"] = counts["source_gene"].map(core_gene)
    counts = counts.dropna(subset=["gene_id"])
    lane_columns = list(range(len(stage_values)))
    for column in lane_columns:
        counts[column] = pd.to_numeric(counts[column], errors="raise")

    lane_groups: dict[str, list[int]] = {}
    records: list[dict[str, object]] = []
    for index, (day, biological_id, lane_id) in enumerate(zip(stage_values, biological_ids, lane_ids, strict=True)):
        sample_id = f"F{day}_{biological_id}"
        lane_groups.setdefault(sample_id, []).append(index)
        records.append({"sample_id": sample_id, "day_after_anthesis": day, "biological_id": biological_id, "lane": lane_id})
    sample_metadata = pd.DataFrame(records).drop_duplicates("sample_id")
    sample_metadata["study_id"] = "GSE128739"
    sample_metadata["genotype"] = "Moneymaker"
    sample_metadata["stage"] = sample_metadata["day_after_anthesis"].map(lambda day: f"day_{day}")
    sample_metadata["stage_ordinal"] = sample_metadata["day_after_anthesis"].map(STAGE_ORDINAL)
    sample_metadata["replicate"] = sample_metadata.groupby("day_after_anthesis").cumcount() + 1
    sample_metadata["include_primary"] = sample_metadata["day_after_anthesis"].isin([41, 53])
    sample_metadata["lane_count"] = sample_metadata["sample_id"].map(lambda sample: len(lane_groups[sample]))
    if (sample_metadata["lane_count"] != 2).any():
        raise ValueError("Expected exactly two sequencing lanes per biological library.")

    collapsed = pd.DataFrame({"gene_id": counts["gene_id"]})
    for sample_id, columns in lane_groups.items():
        collapsed[sample_id] = counts[columns].sum(axis=1)
    collapsed = collapsed.groupby("gene_id", as_index=False).sum()
    libraries = collapsed.drop(columns="gene_id").sum(axis=0)
    if (libraries <= 0).any():
        raise ValueError("A lane-collapsed biological library has non-positive depth.")
    normalized = np.log2(collapsed.drop(columns="gene_id").div(libraries, axis=1) * 1_000_000 + 1)
    normalized.insert(0, "gene_id", collapsed["gene_id"])

    args.output_dir.mkdir(parents=True, exist_ok=True)
    collapsed.to_csv(args.output_dir / "GSE128739_lane_collapsed_counts.csv.gz", index=False, compression="gzip")
    normalized.to_csv(args.output_dir / "GSE128739_log2cpm.csv.gz", index=False, compression="gzip")
    sample_metadata.to_csv(args.output_dir / "GSE128739_sample_metadata.csv", index=False)
    print(f"Prepared GSE128739: {len(collapsed):,} Solyc genes and {len(sample_metadata)} biological libraries; lanes were summed, not replicated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
