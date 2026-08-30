#!/usr/bin/env python3
"""Validate a fixed cross-microarray signature in an independent RNA-seq cohort."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--genotype", default="WT")
    parser.add_argument("--early-ordinal", type=int, default=0)
    parser.add_argument("--late-ordinal", type=int, default=1)
    parser.add_argument("--random-draws", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260826)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    signature = pd.read_csv(args.signature).set_index("gene_id")
    matrix = pd.read_csv(args.matrix, compression="gzip").set_index("gene_id")
    metadata = pd.read_csv(args.metadata)
    metadata = metadata[metadata["genotype"] == args.genotype]
    early = metadata.loc[metadata["stage_ordinal"] == args.early_ordinal, "sample_id"].tolist()
    late = metadata.loc[metadata["stage_ordinal"] == args.late_ordinal, "sample_id"].tolist()
    if len(early) < 2 or len(late) < 2:
        raise ValueError("Each requested endpoint must contain at least two biological libraries.")
    effect = matrix[late].mean(axis=1) - matrix[early].mean(axis=1)
    early_labels = sorted(metadata.loc[metadata["stage_ordinal"] == args.early_ordinal, "stage"].unique())
    late_labels = sorted(metadata.loc[metadata["stage_ordinal"] == args.late_ordinal, "stage"].unique())
    common = signature.index.intersection(effect.index).sort_values()
    source_effect = signature.loc[common, "effect_source_a"]
    validation_effect = effect.loc[common]
    observed_direction = float((np.sign(source_effect) == np.sign(validation_effect)).mean())
    observed_correlation = float(source_effect.rank().corr(validation_effect.rank(), method="pearson"))
    # Permute the measured validation effects across the same genes. This keeps
    # both directional marginals and effect-size distribution fixed under null.
    rng = np.random.default_rng(args.seed)
    random_direction = []
    for _ in range(args.random_draws):
        permuted = rng.permutation(validation_effect.to_numpy())
        random_direction.append(float((np.sign(source_effect.to_numpy()) == np.sign(permuted)).mean()))
    random_array = np.asarray(random_direction)
    report = {
        "scope": "Independent RNA-seq direction check of a signature frozen from two microarray studies; not causal validation.",
        "signature_gene_count_measured": int(len(common)),
        "comparison": f"Imported study endpoint: {', '.join(late_labels)} minus {', '.join(early_labels)}",
        "signature_reference": "GSE42783 red_ripe minus mature_green",
        "spearman_effect_correlation": observed_correlation,
        "direction_concordance": observed_direction,
        "random_direction_concordance": {
            "median": float(np.median(random_array)),
            "upper_95_percentile": float(np.percentile(random_array, 95)),
            "empirical_p": float((1 + np.sum(random_array >= observed_direction)) / (1 + args.random_draws)),
        },
        "limitations": [
            "The imported endpoint labels are study-specific and are not assumed equivalent to the source-A stage labels.",
            "A fixed signature was imported without retraining.",
            "The null permutes validation effects within the same measured signature genes; it does not assess pathway-level confounding.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
