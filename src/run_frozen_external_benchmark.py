#!/usr/bin/env python3
"""Compare frozen development-set baselines on a never-used external cohort."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from studyshield import collapse_independence_groups, evaluate_prediction, random_effects_scores, studyshield_scores


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--development", type=Path, required=True)
    parser.add_argument("--external", type=Path, required=True)
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--biger-scores", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    signature = pd.read_csv(args.signature).set_index("gene_id")
    genes = set(signature.index)
    development = collapse_independence_groups(pd.read_csv(args.development, compression="infer"))
    development = development[development.gene_id.isin(genes)].copy()
    external = pd.read_csv(args.external, compression="infer")
    external = external[external.gene_id.isin(genes)].copy()
    if external.independence_group.nunique() != 1:
        raise ValueError("External benchmark requires exactly one held-out independence group")
    held = external.set_index("gene_id")["signed_effect_rank"]
    matrix = development.pivot(index="gene_id", columns="cohort_id", values="signed_effect_rank")
    shield = studyshield_scores(matrix).set_index("gene_id")
    confidence = matrix.apply(lambda row: max((row.dropna() > 0).mean(), (row.dropna() < 0).mean()) if row.notna().any() else np.nan, axis=1)
    vote = np.sign(matrix).mean(axis=1)
    random_effects = random_effects_scores(development)
    biger = pd.read_csv(args.biger_scores).set_index("gene_id")
    candidates = {
        "studyshield": (shield.studyshield_score, shield.direction_consistency),
        "direction_vote": (vote, confidence),
        "random_effects": (random_effects, confidence),
        "frozen_signature": (signature.effect_source_a, pd.Series(1.0, index=signature.index)),
        "biger_abs_rank_signed_wrapper": (biger.biger_signed_score, biger.direction_probability),
    }
    records = []
    for method, (prediction, certainty) in candidates.items():
        records.append({"held_out_independence_group": external.independence_group.iloc[0], "method": method,
                        **evaluate_prediction(prediction, certainty, held, min(200, len(held)))})
    result = pd.DataFrame(records)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    summary = {
        "external_independence_group": external.independence_group.iloc[0],
        "external_gene_count": int(len(held)),
        "best_method_by_spearman": result.loc[result.spearman_effect_rank.idxmax(), "method"],
        "biger_spearman": float(result.loc[result.method == "biger_abs_rank_signed_wrapper", "spearman_effect_rank"].iloc[0]),
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
