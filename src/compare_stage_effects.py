#!/usr/bin/env python3
"""Compare mature-green to ripe effect directions across two independent studies."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def load(path: Path, metadata_path: Path, early: int, late: int) -> pd.Series:
    matrix = pd.read_csv(path, compression="gzip").set_index("gene_id")
    metadata = pd.read_csv(metadata_path)
    if "genotype" in metadata.columns:
        metadata = metadata[metadata["genotype"] == "WT"]
    early_samples = metadata.loc[metadata["stage_ordinal"] == early, "sample_id"].tolist()
    late_samples = metadata.loc[metadata["stage_ordinal"] == late, "sample_id"].tolist()
    if min(len(early_samples), len(late_samples)) < 2:
        raise ValueError("Each comparison requires at least two samples per stage.")
    return matrix[late_samples].mean(axis=1) - matrix[early_samples].mean(axis=1)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-a-matrix", type=Path, required=True)
    parser.add_argument("--study-a-metadata", type=Path, required=True)
    parser.add_argument("--study-b-matrix", type=Path, required=True)
    parser.add_argument("--study-b-metadata", type=Path, required=True)
    parser.add_argument("--top-k", type=int, default=500)
    parser.add_argument("--random-draws", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260826)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    effect_a = load(args.study_a_matrix, args.study_a_metadata, early=0, late=3)
    effect_b = load(args.study_b_matrix, args.study_b_metadata, early=0, late=3)
    common = effect_a.index.intersection(effect_b.index).sort_values()
    effect_a, effect_b = effect_a.loc[common], effect_b.loc[common]
    # Spearman equals the Pearson correlation of ordinal ranks; avoids a SciPy dependency.
    correlation = float(effect_a.rank(method="average").corr(effect_b.rank(method="average"), method="pearson"))
    top = effect_a.abs().nlargest(args.top_k).index
    observed = float((np.sign(effect_a.loc[top]) == np.sign(effect_b.loc[top])).mean())
    rng = np.random.default_rng(args.seed)
    # Draw one shared random set per iteration so gene directions are compared pairwise.
    background = np.asarray([
        (lambda chosen: (np.sign(effect_a.loc[chosen]) == np.sign(effect_b.loc[chosen])).mean())
        (rng.choice(common, size=args.top_k, replace=False)) for _ in range(args.random_draws)
    ])
    report = {
        "scope": "Cross-study effect-direction check in within-sample rank space; not a causal claim.",
        "comparison": "red_ripe minus mature_green, WT samples only",
        "common_gene_count": int(len(common)),
        "spearman_effect_correlation": correlation,
        "top_k_by_absolute_effect_in_study_a": args.top_k,
        "top_k_direction_concordance": observed,
        "random_direction_concordance": {
            "median": float(np.median(background)),
            "upper_95_percentile": float(np.percentile(background, 95)),
            "empirical_p": float((1 + np.sum(background >= observed)) / (1 + args.random_draws)),
        },
        "limitation": "Feature selection uses study A only; the study-B comparison is external but this is still one pair of studies.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
