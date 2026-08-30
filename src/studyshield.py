"""Core StudyShield estimators and benchmark metrics."""

from __future__ import annotations

import numpy as np
import pandas as pd


WEIGHTS = {"direction": 0.40, "trimmed_mean": 0.30, "weakest": 0.20, "coverage": 0.10}
HETEROGENEITY_PENALTY = 0.50


def collapse_independence_groups(effects: pd.DataFrame) -> pd.DataFrame:
    """Collapse multiple strata from one study before fitting or hold-out evaluation.

    Separate cultivar/tissue contrasts remain visible in the public contrast table,
    but a shared publication/lab family contributes exactly one analysis unit.
    """
    required = {"cohort_id", "independence_group", "gene_id", "signed_effect_rank"}
    missing = required.difference(effects.columns)
    if missing:
        raise ValueError(f"Cannot enforce study independence; missing columns: {sorted(missing)}")
    labels = effects.groupby("independence_group")["cohort_id"].agg(lambda x: "+".join(sorted(set(x))))
    frame = effects.copy()
    frame["study_unit_id"] = frame["independence_group"].map(labels)
    aggregation: dict[str, object] = {"signed_effect_rank": "mean"}
    if "raw_effect" in frame:
        aggregation["raw_effect"] = "mean"
    if "hedges_g" in frame:
        aggregation["hedges_g"] = "mean"
    if "hedges_g_variance" in frame:
        aggregation["hedges_g_variance"] = lambda x: float(x.sum()) / (len(x) ** 2)
    for column in ("n_early", "n_late"):
        if column in frame:
            aggregation[column] = "sum"
    collapsed = frame.groupby(["study_unit_id", "independence_group", "gene_id"], as_index=False).agg(aggregation)
    collapsed = collapsed.rename(columns={"study_unit_id": "cohort_id"})
    return collapsed


def _trimmed_mean(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan
    if len(values) >= 5:
        values = np.sort(values)[1:-1]
    return float(values.mean())


def studyshield_scores(effects: pd.DataFrame, total_cohorts: int | None = None) -> pd.DataFrame:
    """Aggregate a gene x cohort signed-effect-rank matrix without pooling samples."""
    total = total_cohorts or effects.shape[1]
    rows = []
    for gene_id, row in effects.iterrows():
        values = row.dropna().to_numpy(dtype=float)
        if not len(values):
            continue
        positive = int((values > 0).sum())
        negative = int((values < 0).sum())
        direction_consistency = max(positive, negative) / len(values)
        trimmed = _trimmed_mean(values)
        consensus_sign = 1.0 if trimmed >= 0 else -1.0
        weakest = float(np.min(np.abs(values)))
        heterogeneity = float(np.clip(np.std(values, ddof=0), 0, 1))
        coverage = len(values) / total
        magnitude = (
            WEIGHTS["direction"] * direction_consistency
            + WEIGHTS["trimmed_mean"] * abs(trimmed)
            + WEIGHTS["weakest"] * weakest
            + WEIGHTS["coverage"] * coverage
        ) * (1 - HETEROGENEITY_PENALTY * heterogeneity)
        rows.append({
            "gene_id": gene_id,
            "studyshield_score": consensus_sign * magnitude,
            "direction_consistency": direction_consistency,
            "trimmed_mean_effect_rank": trimmed,
            "weakest_absolute_effect_rank": weakest,
            "heterogeneity": heterogeneity,
            "cohort_coverage": coverage,
            "cohort_count": len(values),
        })
    result = pd.DataFrame(rows)
    result["studyshield_rank"] = result["studyshield_score"].abs().rank(method="min", ascending=False).astype(int)
    return result.sort_values("studyshield_rank")


def random_effects_scores(group: pd.DataFrame) -> pd.Series:
    """DerSimonian-Laird random-effects estimate per gene."""
    estimates: dict[str, float] = {}
    for gene_id, rows in group.groupby("gene_id", sort=False):
        y = rows["hedges_g"].to_numpy(float)
        variance = np.maximum(rows["hedges_g_variance"].to_numpy(float), 1e-9)
        if len(y) == 1:
            estimates[gene_id] = float(y[0])
            continue
        weight = 1 / variance
        fixed = np.sum(weight * y) / np.sum(weight)
        q = np.sum(weight * (y - fixed) ** 2)
        c = np.sum(weight) - np.sum(weight**2) / np.sum(weight)
        tau2 = max(0.0, (q - (len(y) - 1)) / c) if c > 0 else 0.0
        random_weight = 1 / (variance + tau2)
        estimates[gene_id] = float(np.sum(random_weight * y) / np.sum(random_weight))
    return pd.Series(estimates, name="random_effects")


def expected_calibration_error(probability: np.ndarray, outcome: np.ndarray, bins: int = 5) -> float:
    edges = np.linspace(0, 1, bins + 1)
    error = 0.0
    for left, right in zip(edges[:-1], edges[1:], strict=True):
        selected = (probability >= left) & (probability <= right if right == 1 else probability < right)
        if selected.any():
            error += selected.mean() * abs(probability[selected].mean() - outcome[selected].mean())
    return float(error)


def evaluate_prediction(prediction: pd.Series, confidence: pd.Series, held_effect: pd.Series, top_k: int) -> dict[str, float]:
    common = prediction.dropna().index.intersection(held_effect.dropna().index)
    predicted = prediction.loc[common]
    observed = held_effect.loc[common]
    correlation = predicted.rank().corr(observed.rank(), method="pearson") if len(common) > 2 else np.nan
    direction = float((np.sign(predicted) == np.sign(observed)).mean())
    k = min(top_k, len(common))
    predicted_top = set(predicted.abs().nlargest(k).index)
    observed_top = set(observed.abs().nlargest(k).index)
    jaccard = len(predicted_top & observed_top) / len(predicted_top | observed_top) if k else np.nan
    outcome = (np.sign(predicted) == np.sign(observed)).astype(float).to_numpy()
    probability = confidence.reindex(common).fillna(0.5).clip(0.5, 1).to_numpy()
    return {
        "gene_count": int(len(common)),
        "spearman_effect_rank": float(correlation),
        "direction_concordance": direction,
        "top_k": int(k),
        "top_k_jaccard": float(jaccard),
        "calibration_error": expected_calibration_error(probability, outcome),
    }
