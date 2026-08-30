#!/usr/bin/env python3
"""Calculate empirical StudyShield p-values from exact within-cohort label permutations."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from studyshield import studyshield_scores


def bh(pvalues: pd.Series) -> pd.Series:
    order = pvalues.sort_values().index
    adjusted = pvalues.loc[order].to_numpy() * len(pvalues) / np.arange(1, len(pvalues) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    return pd.Series(np.minimum(adjusted, 1), index=order).reindex(pvalues.index)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--permutations", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    observed = pd.read_csv(args.observed, compression="infer").set_index("gene_id")
    permutations = pd.read_csv(args.permutations, compression="infer")
    null_scores = []
    for permutation_id, rows in permutations.groupby("permutation_id", sort=True):
        if "independence_group" in rows:
            rows = rows.groupby(["independence_group", "gene_id"], as_index=False)["signed_effect_rank"].mean()
            matrix = rows.pivot(index="gene_id", columns="independence_group", values="signed_effect_rank")
        else:
            matrix = rows.pivot(index="gene_id", columns="cohort_id", values="signed_effect_rank")
        scores = studyshield_scores(matrix)[["gene_id", "studyshield_score"]]
        scores["permutation_id"] = permutation_id
        null_scores.append(scores)
    null = pd.concat(null_scores, ignore_index=True)
    draw_count = permutations.permutation_id.nunique()
    null_matrix = null.pivot(index="gene_id", columns="permutation_id", values="studyshield_score").abs()
    common = observed.index.intersection(null_matrix.index)
    observed_scores = observed.loc[common, "studyshield_score"]
    exceedances = null_matrix.loc[common].ge(observed_scores.abs(), axis=0).sum(axis=1)
    available = null_matrix.loc[common].notna().sum(axis=1)
    result = pd.DataFrame({
        "gene_id": common,
        "studyshield_score": observed_scores.to_numpy(),
        "permutation_count": available.to_numpy(),
        "empirical_p": ((1 + exceedances) / (1 + available)).to_numpy(),
    })
    result["bh_fdr"] = bh(result["empirical_p"])
    result = result.sort_values(["bh_fdr", "empirical_p", "studyshield_score"], ascending=[True, True, False])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"Scored {len(result)} genes against {draw_count} exact within-cohort label permutations.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
