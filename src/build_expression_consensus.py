#!/usr/bin/env python3
"""Build a traceable, expression-only consensus table for a frozen signature.

This script assigns no regulatory edges or causal labels. It simply records
whether each already-frozen gene has a matching effect direction in independent
RNA-seq validations, then ranks genes by the weakest and mean standardized
effect magnitude across sources.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def effect(matrix_path: Path, metadata_path: Path, genotype: str, early: int, late: int) -> pd.Series:
    matrix = pd.read_csv(matrix_path, compression="gzip").set_index("gene_id")
    metadata = pd.read_csv(metadata_path)
    metadata = metadata[metadata["genotype"] == genotype]
    early_samples = metadata.loc[metadata["stage_ordinal"] == early, "sample_id"].tolist()
    late_samples = metadata.loc[metadata["stage_ordinal"] == late, "sample_id"].tolist()
    if len(early_samples) < 2 or len(late_samples) < 2:
        raise ValueError(f"Insufficient biological libraries for {matrix_path}: {early}->{late}")
    return matrix[late_samples].mean(axis=1) - matrix[early_samples].mean(axis=1)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--validation-a-matrix", type=Path, required=True)
    parser.add_argument("--validation-a-metadata", type=Path, required=True)
    parser.add_argument("--validation-a-genotype", required=True)
    parser.add_argument("--validation-a-early", type=int, required=True)
    parser.add_argument("--validation-a-late", type=int, required=True)
    parser.add_argument("--validation-b-matrix", type=Path, required=True)
    parser.add_argument("--validation-b-metadata", type=Path, required=True)
    parser.add_argument("--validation-b-genotype", required=True)
    parser.add_argument("--validation-b-early", type=int, required=True)
    parser.add_argument("--validation-b-late", type=int, required=True)
    parser.add_argument("--gene-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    signature = pd.read_csv(args.signature).set_index("gene_id")
    validation_a = effect(args.validation_a_matrix, args.validation_a_metadata, args.validation_a_genotype, args.validation_a_early, args.validation_a_late)
    validation_b = effect(args.validation_b_matrix, args.validation_b_metadata, args.validation_b_genotype, args.validation_b_early, args.validation_b_late)
    table = signature.join(validation_a.rename("effect_GSE267238"), how="left").join(validation_b.rename("effect_GSE128739"), how="left")
    effect_columns = ["effect_source_a", "effect_source_b", "effect_GSE267238", "effect_GSE128739"]
    table["all_four_measured"] = table[effect_columns].notna().all(axis=1)
    direction = table[effect_columns].apply(lambda column: column * table["effect_source_a"] > 0)
    table["matching_direction_count"] = direction.sum(axis=1)
    table["all_four_same_direction"] = table["all_four_measured"] & direction.all(axis=1)
    standardized = table[effect_columns].abs().rank(pct=True)
    table["minimum_absolute_effect_percentile"] = standardized.min(axis=1)
    table["mean_absolute_effect_percentile"] = standardized.mean(axis=1)
    table["evidence_scope"] = "cross_cohort_expression_replication_only"
    table["candidate_rank"] = pd.NA
    accepted = table["all_four_same_direction"]
    table.loc[accepted, "candidate_rank"] = (
        table.loc[accepted]
        .sort_values(["minimum_absolute_effect_percentile", "mean_absolute_effect_percentile"], ascending=False)
        .assign(candidate_rank=lambda frame: range(1, len(frame) + 1))["candidate_rank"]
    )
    annotation = pd.read_csv(args.gene_metadata, compression="gzip").drop_duplicates("gene_id")
    result = table.reset_index().merge(annotation, on="gene_id", how="left")
    result = result.sort_values(["all_four_same_direction", "candidate_rank"], ascending=[False, True])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(
        f"Wrote {len(result)} frozen-signature genes; "
        f"{int(result['all_four_same_direction'].sum())} match direction in all four expression sources."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
