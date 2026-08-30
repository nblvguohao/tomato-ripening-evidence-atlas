#!/usr/bin/env python3
"""Fit StudyShield, bootstrap stability, and leave-one-cohort-out benchmarks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from studyshield import collapse_independence_groups, evaluate_prediction, random_effects_scores, studyshield_scores


def baseline_predictions(train: pd.DataFrame, frozen: pd.DataFrame) -> dict[str, tuple[pd.Series, pd.Series]]:
    rank = train.pivot(index="gene_id", columns="cohort_id", values="signed_effect_rank")
    shield = studyshield_scores(rank).set_index("gene_id")
    vote = np.sign(rank).mean(axis=1)
    weighted = train.assign(weight=train["n_early"] + train["n_late"])
    pooled = weighted.assign(weighted_effect=weighted["signed_effect_rank"] * weighted["weight"]).groupby("gene_id").agg(
        numerator=("weighted_effect", "sum"), denominator=("weight", "sum"), confidence=("signed_effect_rank", lambda x: max((x > 0).mean(), (x < 0).mean()))
    )
    pooled_score = pooled["numerator"] / pooled["denominator"]
    random_effects = random_effects_scores(train)
    random_confidence = rank.apply(lambda row: max((row.dropna() > 0).mean(), (row.dropna() < 0).mean()) if row.notna().any() else np.nan, axis=1)
    frozen_score = frozen.set_index("gene_id")["effect_source_a"]
    frozen_confidence = pd.Series(1.0, index=frozen_score.index)
    return {
        "studyshield": (shield["studyshield_score"], shield["direction_consistency"]),
        "direction_vote": (vote, random_confidence),
        "random_effects": (random_effects, random_confidence),
        "pooled_batch_adjusted": (pooled_score, pooled["confidence"]),
        "frozen_signature": (frozen_score, frozen_confidence),
    }


def bootstrap(effects: pd.DataFrame, draws: int, top_k: int, seed: int) -> pd.DataFrame:
    matrix = effects.pivot(index="gene_id", columns="cohort_id", values="signed_effect_rank")
    cohorts = matrix.columns.to_numpy()
    rng = np.random.default_rng(seed)
    records = []
    for draw in range(draws):
        sampled = rng.choice(cohorts, size=len(cohorts), replace=True)
        sampled_matrix = matrix.loc[:, sampled].copy()
        sampled_matrix.columns = [f"sample_{index}" for index in range(len(sampled))]
        score = studyshield_scores(sampled_matrix).set_index("gene_id")["studyshield_score"]
        rank = score.abs().rank(method="min", ascending=False)
        records.append(pd.DataFrame({"gene_id": score.index, "draw": draw + 1, "score": score.values, "rank": rank.values, "selected_top_k": rank.values <= top_k}))
    all_draws = pd.concat(records, ignore_index=True)
    return all_draws.groupby("gene_id", as_index=False).agg(
        bootstrap_score_median=("score", "median"),
        bootstrap_score_lower_95=("score", lambda x: x.quantile(.025)),
        bootstrap_score_upper_95=("score", lambda x: x.quantile(.975)),
        bootstrap_rank_median=("rank", "median"),
        top_k_selection_probability=("selected_top_k", "mean"),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contrasts", type=Path, required=True)
    parser.add_argument("--frozen-signature", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap-draws", type=int, default=500)
    parser.add_argument("--top-k", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260827)
    args = parser.parse_args()
    contrast_strata = pd.read_csv(args.contrasts, compression="infer")
    effects = collapse_independence_groups(contrast_strata)
    frozen = pd.read_csv(args.frozen_signature)
    signature_genes = set(frozen["gene_id"])
    evaluation = effects[effects["gene_id"].isin(signature_genes)].copy()
    matrix = effects.pivot(index="gene_id", columns="cohort_id", values="signed_effect_rank")
    fitted = studyshield_scores(matrix)
    stability = bootstrap(effects[effects["gene_id"].isin(signature_genes)], args.bootstrap_draws, min(args.top_k, len(signature_genes)), args.seed)
    fitted = fitted.merge(stability, on="gene_id", how="left", validate="one_to_one")

    benchmark_rows = []
    for held_cohort in sorted(evaluation["cohort_id"].unique()):
        train = evaluation[evaluation["cohort_id"] != held_cohort]
        training_cohorts = sorted(train["cohort_id"].unique())
        if held_cohort in training_cohorts:
            raise AssertionError(f"LOCO leakage detected for {held_cohort}")
        held = evaluation[evaluation["cohort_id"] == held_cohort].set_index("gene_id")["signed_effect_rank"]
        for method, (prediction, confidence) in baseline_predictions(train, frozen).items():
            metrics = evaluate_prediction(prediction, confidence, held, min(args.top_k, len(held)))
            benchmark_rows.append({
                "held_out_cohort": held_cohort,
                "training_cohort_count": len(training_cohorts),
                "training_cohorts": ";".join(training_cohorts),
                "method": method,
                **metrics,
            })
    benchmark = pd.DataFrame(benchmark_rows)
    competitors = ["direction_vote", "random_effects", "pooled_batch_adjusted", "frozen_signature"]
    wins = []
    for cohort, rows in benchmark.groupby("held_out_cohort"):
        shield = rows[rows.method == "studyshield"].iloc[0]
        wins.append(all(shield["spearman_effect_rank"] > rows[rows.method == method]["spearman_effect_rank"].iloc[0] for method in competitors))
    summary = {
        "scope": "Study-level cross-platform aggregation; samples are never pooled across cohorts.",
        "contrast_strata_count": int(contrast_strata.cohort_id.nunique()),
        "cohort_count": int(effects.cohort_id.nunique()),
        "independence_group_count": int(effects.independence_group.nunique()),
        "fitted_gene_count": int(len(fitted)),
        "frozen_signature_gene_count": int(len(signature_genes)),
        "bootstrap_draws": args.bootstrap_draws,
        "top_k": args.top_k,
        "studyshield_beats_all_baselines_by_spearman_in_cohorts": int(sum(wins)),
        "held_out_cohort_count": len(wins),
        "method_superiority_gate_passed": bool(sum(wins) > len(wins) / 2),
        "claim_boundary": "If the superiority gate fails, StudyShield is retained only as an atlas integration score.",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fitted.to_csv(args.output_dir / "studyshield_scores.csv.gz", index=False, compression="gzip")
    benchmark.to_csv(args.output_dir / "studyshield_loco_benchmark.csv", index=False)
    (args.output_dir / "studyshield_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
