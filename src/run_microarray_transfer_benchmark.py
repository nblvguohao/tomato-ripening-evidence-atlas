#!/usr/bin/env python3
"""Frozen external stage-transfer benchmark in shared, within-sample rank space.

The method is intentionally simple: ANOVA feature selection is fitted only on
the training study, followed by a nearest-centroid classifier. Random-feature
controls show whether selected features outperform equally sized random sets.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


LABELS = {0: "mature_green", 1: "transition", 3: "red_ripe"}


def load_matrix(path: Path, metadata: pd.DataFrame) -> pd.DataFrame:
    matrix = pd.read_csv(path, compression="gzip").set_index("gene_id")
    return matrix.loc[:, metadata["sample_id"].tolist()]


def stage_labels(metadata: pd.DataFrame, training: bool) -> np.ndarray:
    values = metadata["stage_ordinal"].to_numpy()
    if training:
        return values
    return values  # GSE108415 uses ordinal 0, 1, 3 by construction


def anova_f_values(features: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """One-way ANOVA F statistic calculated without a third-party ML package."""
    classes = np.unique(labels)
    overall = features.mean(axis=0)
    between = np.zeros(features.shape[1])
    within = np.zeros(features.shape[1])
    for label in classes:
        subset = features[labels == label]
        between += len(subset) * (subset.mean(axis=0) - overall) ** 2
        within += ((subset - subset.mean(axis=0)) ** 2).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (between / (len(classes) - 1)) / (within / (len(features) - len(classes)))


def predict_nearest_centroid(train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray) -> np.ndarray:
    classes = np.unique(train_y)
    centroids = np.vstack([train_x[train_y == label].mean(axis=0) for label in classes])
    squared_distance = ((test_x[:, None, :] - centroids[None, :, :]) ** 2).sum(axis=2)
    return classes[squared_distance.argmin(axis=1)]


def confusion_and_metrics(actual: np.ndarray, prediction: np.ndarray) -> tuple[list[list[int]], float, float]:
    labels = [0, 1, 3]
    matrix = np.array([[int(np.sum((actual == a) & (prediction == p))) for p in labels] for a in labels])
    per_class_f1 = []
    for index in range(len(labels)):
        tp = matrix[index, index]
        fp = matrix[:, index].sum() - tp
        fn = matrix[index, :].sum() - tp
        denominator = 2 * tp + fp + fn
        per_class_f1.append((2 * tp / denominator) if denominator else 0.0)
    return matrix.tolist(), float(np.mean(actual == prediction)), float(np.mean(per_class_f1))


def score_features(train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray, test_y: np.ndarray, indices: np.ndarray) -> dict[str, object]:
    prediction = predict_nearest_centroid(train_x[:, indices], train_y, test_x[:, indices])
    matrix, accuracy, macro_f1 = confusion_and_metrics(test_y, prediction)
    return {
        "accuracy": accuracy,
        "macro_f1": macro_f1,
        "confusion_matrix": matrix,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--train-matrix", type=Path, required=True)
    parser.add_argument("--train-metadata", type=Path, required=True)
    parser.add_argument("--test-matrix", type=Path, required=True)
    parser.add_argument("--test-metadata", type=Path, required=True)
    parser.add_argument("--top-k", type=int, default=150)
    parser.add_argument("--random-draws", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260826)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    train_meta = pd.read_csv(args.train_metadata)
    test_meta = pd.read_csv(args.test_metadata)
    test_meta = test_meta[test_meta["genotype"] == "WT"].copy()  # preserve external genotype shift as a later sensitivity test
    train_matrix = load_matrix(args.train_matrix, train_meta)
    test_matrix = load_matrix(args.test_matrix, test_meta)
    common_genes = train_matrix.index.intersection(test_matrix.index).sort_values()
    if len(common_genes) < args.top_k:
        raise ValueError("Too few common genes for requested top-k.")

    train_x = train_matrix.loc[common_genes].T.to_numpy()
    test_x = test_matrix.loc[common_genes].T.to_numpy()
    train_y = stage_labels(train_meta, training=True)
    test_y = stage_labels(test_meta, training=False)
    f_values = anova_f_values(train_x, train_y)
    ranked_indices = np.argsort(np.nan_to_num(f_values, nan=-np.inf))[::-1]
    selected = ranked_indices[:args.top_k]
    main_result = score_features(train_x, train_y, test_x, test_y, selected)

    rng = np.random.default_rng(args.seed)
    random_f1 = [score_features(train_x, train_y, test_x, test_y, rng.choice(len(common_genes), size=args.top_k, replace=False))["macro_f1"]
                 for _ in range(args.random_draws)]
    random_f1_array = np.asarray(random_f1, dtype=float)
    report = {
        "scope": "Probe/gene-label transfer benchmark only; not a causal-network result.",
        "train_study": "GSE42783",
        "test_study": "GSE108415_WT_only",
        "common_gene_count": int(len(common_genes)),
        "top_k": args.top_k,
        "stage_labels": LABELS,
        "selected_feature_result": main_result,
        "random_feature_macro_f1": {
            "draw_count": args.random_draws,
            "median": float(np.median(random_f1_array)),
            "upper_95_percentile": float(np.percentile(random_f1_array, 95)),
            "empirical_p": float((1 + np.sum(random_f1_array >= main_result["macro_f1"])) / (1 + args.random_draws)),
        },
        "top_features": common_genes[selected].tolist(),
        "limitations": [
            "Training study has only three samples per class.",
            "Feature count was prespecified, not tuned on the external study.",
            "Within-sample rank transformation removes absolute expression scale.",
            "The test is WT-only; LCY transgenic samples are reserved for a separate perturbation sensitivity analysis.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
