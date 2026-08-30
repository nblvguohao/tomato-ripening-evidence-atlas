#!/usr/bin/env python3
"""Create one standardized gene-effect table across independently processed cohorts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_CONFIG = {
    "cohort_id", "matrix_path", "metadata_path", "assay", "cultivar", "tissue",
    "filter_column", "filter_value", "early_ordinal", "late_ordinal",
    "independence_group", "analysis_role", "admission_status",
}


def hedges_g_and_variance(early: np.ndarray, late: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n_early, n_late = early.shape[1], late.shape[1]
    mean_difference = late.mean(axis=1) - early.mean(axis=1)
    pooled_variance = ((n_early - 1) * early.var(axis=1, ddof=1) + (n_late - 1) * late.var(axis=1, ddof=1)) / (n_early + n_late - 2)
    pooled_sd = np.sqrt(pooled_variance)
    raw_d = np.divide(mean_difference, pooled_sd, out=np.zeros_like(mean_difference), where=pooled_sd > 0)
    correction = 1 - 3 / (4 * (n_early + n_late) - 9)
    g = correction * raw_d
    variance = (n_early + n_late) / (n_early * n_late) + np.divide(g * g, 2 * (n_early + n_late - 2), out=np.zeros_like(g), where=(n_early + n_late - 2) > 0)
    return g, variance


def signed_absolute_rank(effect: pd.Series) -> pd.Series:
    magnitude = effect.abs().rank(method="average", pct=True)
    return np.sign(effect) * magnitude


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--signature", type=Path)
    parser.add_argument("--permutation-output", type=Path)
    parser.add_argument("--permutations", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260827)
    args = parser.parse_args()
    if bool(args.permutation_output) != bool(args.signature and args.permutations > 0):
        raise ValueError("Permutation output requires --signature and a positive --permutations value")
    signature_genes = set(pd.read_csv(args.signature)["gene_id"]) if args.signature else set()
    rng = np.random.default_rng(args.seed)
    config = pd.read_csv(args.config, dtype=str).fillna("")
    missing = REQUIRED_CONFIG.difference(config.columns)
    if missing:
        raise ValueError(f"Contrast config is missing columns: {sorted(missing)}")
    if config["cohort_id"].duplicated().any():
        raise ValueError("cohort_id must be unique")

    tables: list[pd.DataFrame] = []
    audits: list[dict[str, object]] = []
    permutation_tables: list[pd.DataFrame] = []
    for row in config.to_dict("records"):
        if row["admission_status"] != "eligible":
            continue
        matrix_path = args.root / row["matrix_path"]
        metadata_path = args.root / row["metadata_path"]
        matrix = pd.read_csv(matrix_path, compression="infer").set_index("gene_id")
        metadata = pd.read_csv(metadata_path)
        if row["filter_column"] not in metadata:
            raise ValueError(f"{row['cohort_id']}: missing filter column {row['filter_column']}")
        selected = metadata[metadata[row["filter_column"]].astype(str) == row["filter_value"]].copy()
        early_ordinal, late_ordinal = int(row["early_ordinal"]), int(row["late_ordinal"])
        early_ids = selected.loc[selected["stage_ordinal"] == early_ordinal, "sample_id"].tolist()
        late_ids = selected.loc[selected["stage_ordinal"] == late_ordinal, "sample_id"].tolist()
        if len(early_ids) < 2 or len(late_ids) < 2:
            raise ValueError(f"{row['cohort_id']}: each endpoint requires at least two biological replicates")
        requested = early_ids + late_ids
        missing_samples = set(requested).difference(matrix.columns)
        if missing_samples:
            raise ValueError(f"{row['cohort_id']}: matrix lacks samples {sorted(missing_samples)}")
        values = matrix[requested].apply(pd.to_numeric, errors="raise")
        effect = values[late_ids].mean(axis=1) - values[early_ids].mean(axis=1)
        g, g_variance = hedges_g_and_variance(values[early_ids].to_numpy(), values[late_ids].to_numpy())
        table = pd.DataFrame({
            "cohort_id": row["cohort_id"],
            "gene_id": values.index,
            "raw_effect": effect.to_numpy(),
            "signed_effect_rank": signed_absolute_rank(effect).to_numpy(),
            "hedges_g": g,
            "hedges_g_variance": g_variance,
            "measured": True,
            "n_early": len(early_ids),
            "n_late": len(late_ids),
            "assay": row["assay"],
            "cultivar": row["cultivar"],
            "tissue": row["tissue"],
            "early_ordinal": early_ordinal,
            "late_ordinal": late_ordinal,
            "independence_group": row["independence_group"],
            "analysis_role": row["analysis_role"],
        })
        tables.append(table)
        if args.permutation_output:
            measured_signature = values.index.intersection(sorted(signature_genes))
            for permutation_id in range(1, args.permutations + 1):
                permuted = rng.permutation(requested)
                permuted_early = list(permuted[:len(early_ids)])
                permuted_late = list(permuted[len(early_ids):])
                permuted_effect = values[permuted_late].mean(axis=1) - values[permuted_early].mean(axis=1)
                permuted_rank = signed_absolute_rank(permuted_effect).reindex(measured_signature)
                permutation_tables.append(pd.DataFrame({
                    "permutation_id": permutation_id,
                    "cohort_id": row["cohort_id"],
                    "independence_group": row["independence_group"],
                    "gene_id": measured_signature,
                    "signed_effect_rank": permuted_rank.to_numpy(),
                }))
        audits.append({
            "cohort_id": row["cohort_id"], "genes": len(table), "early_samples": early_ids,
            "late_samples": late_ids, "assay": row["assay"], "cultivar": row["cultivar"],
            "tissue": row["tissue"], "matrix_path": str(matrix_path), "metadata_path": str(metadata_path),
        })
    output = pd.concat(tables, ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False, compression="gzip" if args.output.suffix == ".gz" else None)
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    args.audit.write_text(json.dumps({"cohort_count": len(audits), "cohorts": audits}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if args.permutation_output:
        permutations = pd.concat(permutation_tables, ignore_index=True)
        args.permutation_output.parent.mkdir(parents=True, exist_ok=True)
        permutations.to_csv(args.permutation_output, index=False, compression="gzip" if args.permutation_output.suffix == ".gz" else None)
        print(f"Wrote {len(permutations):,} exact within-cohort label-permutation effects.")
    print(f"Wrote {len(output):,} cohort-gene effects from {len(audits)} eligible cohorts.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
