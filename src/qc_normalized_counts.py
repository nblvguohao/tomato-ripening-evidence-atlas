#!/usr/bin/env python3
"""Generate transparent sample-level QC for normalized count matrices."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--group-columns", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, compression="gzip").set_index("gene_id")
    metadata = pd.read_csv(args.metadata)
    sample_ids = metadata["sample_id"].tolist()
    if set(sample_ids) != set(matrix.columns):
        raise ValueError("Matrix columns and metadata sample_id values disagree.")

    matrix = matrix.loc[:, sample_ids]
    sample_medians = matrix.median(axis=0)
    correlations = matrix.corr(method="spearman")
    groups = metadata.groupby(args.group_columns, dropna=False, sort=True)["sample_id"].apply(list)
    within_group: list[dict[str, object]] = []
    for labels, members in groups.items():
        label_values = labels if isinstance(labels, tuple) else (labels,)
        pair_values = [float(correlations.loc[a, b]) for a, b in itertools.combinations(members, 2)]
        row = dict(zip(args.group_columns, label_values))
        row.update({"n_samples": len(members), "n_pairs": len(pair_values)})
        row["median_spearman"] = float(pd.Series(pair_values).median()) if pair_values else None
        row["minimum_spearman"] = float(pd.Series(pair_values).min()) if pair_values else None
        within_group.append(row)

    report = {
        "matrix": str(args.matrix),
        "gene_count": int(matrix.shape[0]),
        "sample_count": int(matrix.shape[1]),
        "per_sample_log2cpm_median": {sample: float(value) for sample, value in sample_medians.items()},
        "within_group_spearman": within_group,
        "interpretation": "High within-group correlation supports technical consistency but does not establish biological validity or remove batch effects.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
