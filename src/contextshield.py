"""Context-aware robust aggregation for signed cross-study gene effects.

This is intentionally an empirical-Bayes approximation rather than a claim of
having reproduced BiGER.  It consumes within-study signed effect ranks, never
pooled expression matrices.  Study reliability, context offsets, and
observation-level Huber weights are fitted using only the supplied studies.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ContextShieldConfig:
    iterations: int = 6
    context_shrinkage: float = 3.0
    huber_k: float = 1.75
    min_reliability: float = 0.05


def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not valid.any():
        return np.nan
    return float(np.average(values[valid], weights=weights[valid]))


def _study_quality(frame: pd.DataFrame) -> pd.Series:
    """Return a bounded, evidence-quality prior for each study unit."""
    sizes = frame.groupby("cohort_id")[["n_early", "n_late"]].first().min(axis=1).clip(lower=1)
    coverage = frame.groupby("cohort_id")["measured"].mean().reindex(sizes.index).fillna(0.0)
    quality = np.sqrt(sizes.to_numpy(float)) * coverage.to_numpy(float)
    quality = quality / np.nanmedian(quality) if np.nanmedian(quality) > 0 else np.ones_like(quality)
    return pd.Series(np.clip(quality, 0.25, 4.0), index=sizes.index, name="quality_prior")


def collapse_context_independence_groups(effects: pd.DataFrame) -> pd.DataFrame:
    """Collapse shared-publication strata while retaining reproducible context."""
    required = {"cohort_id", "independence_group", "gene_id", "signed_effect_rank"}
    missing = required.difference(effects.columns)
    if missing:
        raise ValueError(f"Missing required ContextShield columns: {sorted(missing)}")
    copy = effects.copy()
    labels = copy.groupby("independence_group")["cohort_id"].agg(lambda x: "+".join(sorted(set(x))))
    copy["cohort_id"] = copy["independence_group"].map(labels)
    aggregations: dict[str, object] = {"signed_effect_rank": "mean"}
    for numeric in ("raw_effect", "hedges_g", "hedges_g_variance", "n_early", "n_late", "measured"):
        if numeric in copy:
            aggregations[numeric] = "mean" if numeric not in {"n_early", "n_late"} else "sum"
    collapsed = copy.groupby(["cohort_id", "independence_group", "gene_id"], as_index=False).agg(aggregations)
    metadata = effects.groupby("independence_group", as_index=False).agg(
        tissue=("tissue", lambda x: "mixed" if x.nunique() > 1 else str(x.iloc[0])),
        assay=("assay", lambda x: "mixed" if x.nunique() > 1 else str(x.iloc[0])),
        early_ordinal=("early_ordinal", lambda x: "mixed" if x.nunique() > 1 else str(x.iloc[0])),
        late_ordinal=("late_ordinal", lambda x: "mixed" if x.nunique() > 1 else str(x.iloc[0])),
    )
    metadata["transition"] = metadata["early_ordinal"] + "->" + metadata["late_ordinal"]
    return collapsed.merge(metadata.drop(columns=["early_ordinal", "late_ordinal"]), on="independence_group", how="left", validate="many_to_one")


def fit_contextshield(
    effects: pd.DataFrame,
    config: ContextShieldConfig = ContextShieldConfig(),
    context_column: str = "transition",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit core signed effects and return gene and study diagnostics.

    ``effects`` must contain one independent study unit per cohort_id.  All
    reliability estimates are internal to this frame, making it safe for a
    leave-one-independent-group-out caller to avoid leakage.
    """
    required = {"cohort_id", "gene_id", "signed_effect_rank", context_column}
    missing = required.difference(effects.columns)
    if missing:
        raise ValueError(f"Missing required ContextShield columns: {sorted(missing)}")
    frame = effects.copy()
    frame["signed_effect_rank"] = pd.to_numeric(frame["signed_effect_rank"], errors="coerce")
    frame = frame.dropna(subset=["signed_effect_rank"])
    if frame.empty:
        raise ValueError("No finite signed effect ranks available")
    if "measured" not in frame:
        frame["measured"] = True
    if "n_early" not in frame:
        frame["n_early"] = 2
    if "n_late" not in frame:
        frame["n_late"] = 2
    frame[context_column] = frame[context_column].fillna("unknown").astype(str)

    studies = pd.Index(sorted(frame["cohort_id"].unique()), name="cohort_id")
    genes = pd.Index(sorted(frame["gene_id"].unique()), name="gene_id")
    contexts = pd.Index(sorted(frame[context_column].unique()), name=context_column)
    g_index = {gene: i for i, gene in enumerate(genes)}
    s_index = {study: i for i, study in enumerate(studies)}
    c_index = {context: i for i, context in enumerate(contexts)}
    values = np.full((len(genes), len(studies)), np.nan)
    study_context = np.zeros(len(studies), dtype=int)
    for row in frame[["gene_id", "cohort_id", "signed_effect_rank", context_column]].itertuples(index=False):
        values[g_index[row.gene_id], s_index[row.cohort_id]] = row.signed_effect_rank
        study_context[s_index[row.cohort_id]] = c_index[getattr(row, context_column)]

    quality = _study_quality(frame).reindex(studies).fillna(1.0).to_numpy(float)
    theta = np.nanmedian(values, axis=1)
    theta = np.nan_to_num(theta, nan=0.0)
    offsets = np.zeros((len(genes), len(contexts)))
    reliability = quality / quality.sum()
    observation_weight = np.where(np.isfinite(values), 1.0, 0.0)

    for _ in range(config.iterations):
        predicted = theta[:, None] + offsets[:, study_context]
        residual = values - predicted
        study_scale = np.nanmedian(np.abs(residual), axis=0) / 0.6745
        fallback = np.nanmedian(study_scale[np.isfinite(study_scale) & (study_scale > 1e-6)])
        fallback = float(fallback) if np.isfinite(fallback) else 0.25
        study_scale = np.where(np.isfinite(study_scale) & (study_scale > 1e-6), study_scale, fallback)
        standardized = np.abs(residual) / study_scale[None, :]
        observation_weight = np.where(
            np.isfinite(values), np.minimum(1.0, config.huber_k / np.maximum(standardized, 1e-12)), 0.0
        )
        outlier_fraction = np.nanmean(standardized > 2.5, axis=0)

        agreement = []
        for index in range(len(studies)):
            others = np.delete(values, index, axis=1)
            reference = np.nanmedian(others, axis=1)
            observed = values[:, index]
            valid = np.isfinite(reference) & np.isfinite(observed)
            corr = pd.Series(observed[valid]).corr(pd.Series(reference[valid]), method="spearman") if valid.sum() >= 10 else 0.0
            agreement.append(max(0.02, float(corr) if np.isfinite(corr) else 0.02))
        reliability_raw = quality * np.asarray(agreement) * np.clip(1.0 - outlier_fraction, 0.05, 1.0)
        reliability_raw = np.maximum(reliability_raw, config.min_reliability)
        reliability = reliability_raw / reliability_raw.sum()

        combined = observation_weight * reliability[None, :]
        adjusted = values - offsets[:, study_context]
        denominator = np.where(np.isfinite(adjusted), combined, 0.0).sum(axis=1)
        theta = np.divide(
            np.nansum(adjusted * combined, axis=1), denominator,
            out=np.zeros(len(genes)), where=denominator > 0,
        )
        for context_index in range(len(contexts)):
            selected = study_context == context_index
            residual_context = values[:, selected] - theta[:, None]
            context_weight = combined[:, selected]
            n_effective = np.where(np.isfinite(residual_context), context_weight, 0.0).sum(axis=1)
            raw_offset = np.divide(
                np.nansum(residual_context * context_weight, axis=1), n_effective,
                out=np.zeros(len(genes)), where=n_effective > 0,
            )
            offsets[:, context_index] = raw_offset * n_effective / (n_effective + config.context_shrinkage)

    predicted = theta[:, None] + offsets[:, study_context]
    residual = values - predicted
    observed_weight = observation_weight * reliability[None, :]
    positive_weight = (observed_weight * (values > 0)).sum(axis=1)
    negative_weight = (observed_weight * (values < 0)).sum(axis=1)
    total_weight = positive_weight + negative_weight
    direction_probability = np.divide(np.maximum(positive_weight, negative_weight), total_weight, out=np.full(len(genes), 0.5), where=total_weight > 0)
    coverage = np.isfinite(values).mean(axis=1)
    heterogeneity = np.nanstd(residual, axis=1)
    score = theta * (0.5 + 0.5 * direction_probability) * coverage
    genes_result = pd.DataFrame({
        "gene_id": genes,
        "contextshield_score": score,
        "core_signed_effect_rank": theta,
        "direction_probability": direction_probability,
        "cohort_coverage": coverage,
        "residual_heterogeneity": heterogeneity,
        "cohort_count": np.isfinite(values).sum(axis=1),
    })
    for context_index, context in enumerate(contexts):
        genes_result[f"context_offset__{context}"] = offsets[:, context_index]
    genes_result["contextshield_rank"] = genes_result["contextshield_score"].abs().rank(method="min", ascending=False).astype(int)
    genes_result = genes_result.sort_values("contextshield_rank").reset_index(drop=True)
    study_result = pd.DataFrame({
        "cohort_id": studies,
        "context": contexts.take(study_context),
        "quality_prior": quality,
        "reliability": reliability,
        "agreement_spearman": np.asarray(agreement),
        "outlier_fraction": outlier_fraction,
        "residual_scale": study_scale,
    }).sort_values("reliability", ascending=False).reset_index(drop=True)
    return genes_result, study_result
