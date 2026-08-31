#!/usr/bin/env python3
"""Evaluate frozen cross-layer enrichments against effect-matched null genes.

The frozen signature was selected using GSE42783 effects and GSE108415 direction
confirmation. This historical sensitivity script instead matches GSE267238 WT
stage effect, endpoint detection and gene length. It does NOT control the source
selection variables. Selection-support diagnostics are provided separately in
the post-review revision. No signature or assay admission rule is changed here.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


SOLYC = re.compile(r"(Solyc\d{2}g\d{6})")


def core(value: object) -> str:
    match = SOLYC.search(str(value))
    return match.group(1) if match else ""


def truth(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def logcomb(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def conditional_pmf(a: int, b: int, c: int, d: int, log_odds: float) -> tuple[np.ndarray, np.ndarray]:
    """Fisher's noncentral hypergeometric distribution with fixed margins."""
    row1, col1, total = a + b, a + c, a + b + c + d
    lower, upper = max(0, row1 + col1 - total), min(row1, col1)
    values = np.arange(lower, upper + 1)
    weights = np.array([
        logcomb(col1, int(x)) + logcomb(total - col1, row1 - int(x)) + int(x) * log_odds
        for x in values
    ])
    weights = np.exp(weights - weights.max())
    return values, weights / weights.sum()


def exact_conditional_odds_ratio_ci(a: int, b: int, c: int, d: int, alpha: float = 0.05) -> tuple[float, float]:
    """Two-sided exact conditional CI by inverting Fisher's noncentral law."""
    row1, col1, total = a + b, a + c, a + b + c + d
    lower_support, upper_support = max(0, row1 + col1 - total), min(row1, col1)
    if a == lower_support:
        lower = 0.0
    else:
        lo, hi = -50.0, 50.0
        for _ in range(100):
            mid = (lo + hi) / 2
            values, probabilities = conditional_pmf(a, b, c, d, mid)
            if probabilities[values >= a].sum() > alpha / 2:
                hi = mid
            else:
                lo = mid
        lower = math.exp((lo + hi) / 2)
    if a == upper_support:
        upper = math.inf
    else:
        lo, hi = -50.0, 50.0
        for _ in range(100):
            mid = (lo + hi) / 2
            values, probabilities = conditional_pmf(a, b, c, d, mid)
            if probabilities[values <= a].sum() > alpha / 2:
                lo = mid
            else:
                hi = mid
        upper = math.exp((lo + hi) / 2)
    return lower, upper


def fisher_greater(a: int, b: int, c: int, d: int) -> float:
    total, row1, col1 = a + b + c + d, a + b, a + c
    denominator = math.comb(total, row1)
    return float(sum(
        math.comb(col1, x) * math.comb(total - col1, row1 - x) / denominator
        for x in range(a, min(row1, col1) + 1)
    ))


def raw_table(signature: set[str], universe: set[str], target: set[str]) -> dict[str, object]:
    sig = signature & universe
    target = target & universe
    a, b = len(sig & target), len(sig - target)
    c, d = len(target - sig), len(universe - sig - target)
    odds_ratio = (a * d / (b * c)) if b and c else math.inf
    lower, upper = exact_conditional_odds_ratio_ci(a, b, c, d)
    return {
        "signature_target": a, "signature_not_target": b,
        "nonsignature_target": c, "nonsignature_not_target": d,
        "odds_ratio": odds_ratio,
        "exact_conditional_95_ci": [lower, upper],
        "fisher_greater_p": fisher_greater(a, b, c, d),
    }


def gene_lengths(path: Path) -> pd.Series:
    lengths: dict[str, int] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            gene = core(fields[8])
            if gene:
                lengths[gene] = abs(int(fields[4]) - int(fields[3])) + 1
    return pd.Series(lengths, name="gene_length_bp")


def matching_covariates(matrix_path: Path, metadata_path: Path, gff3: Path) -> pd.DataFrame:
    matrix = pd.read_csv(matrix_path, compression="infer").set_index("gene_id")
    matrix.index = matrix.index.map(core)
    matrix = matrix[matrix.index != ""].groupby(level=0).mean()
    metadata = pd.read_csv(metadata_path)
    wt = metadata[metadata["genotype"].eq("WT")]
    early = wt.loc[wt["stage_ordinal"].eq(0), "sample_id"].tolist()
    late = wt.loc[wt["stage_ordinal"].eq(1), "sample_id"].tolist()
    if len(early) < 2 or len(late) < 2:
        raise ValueError("GSE267238 WT endpoints need at least two biological libraries each")
    covariates = pd.DataFrame(index=matrix.index)
    covariates["detection_rate"] = (matrix[early + late] > 1).mean(axis=1)
    covariates["absolute_effect"] = (matrix[late].mean(axis=1) - matrix[early].mean(axis=1)).abs()
    covariates = covariates.join(gene_lengths(gff3), how="inner")
    return covariates[covariates["gene_length_bp"].gt(0)].copy()


def match_layer(
    name: str, signature: set[str], universe: set[str], target: set[str], covariates: pd.DataFrame,
    draws: int, seed: int,
) -> tuple[dict[str, object], pd.DataFrame]:
    eligible_signature = sorted(signature & universe & set(covariates.index))
    candidate_ids = np.array(sorted((universe & set(covariates.index)) - set(eligible_signature)))
    if len(eligible_signature) < 5 or len(candidate_ids) < len(eligible_signature):
        raise ValueError(f"{name}: insufficient covariate-complete assay background for matching")
    values = covariates.loc[list(set(eligible_signature) | set(candidate_ids)), ["detection_rate", "absolute_effect", "gene_length_bp"]].copy()
    values["gene_length_bp"] = np.log10(values["gene_length_bp"])
    z = (values - values.mean()) / values.std().replace(0, 1)
    candidate_z = z.loc[candidate_ids].to_numpy()
    selected_z = z.loc[eligible_signature].to_numpy()
    pool_size = min(200, len(candidate_ids))
    neighbor_pools: list[np.ndarray] = []
    for vector in selected_z:
        distances = ((candidate_z - vector) ** 2).sum(axis=1)
        pool = np.argpartition(distances, pool_size - 1)[:pool_size]
        neighbor_pools.append(pool[np.argsort(distances[pool])])
    rng = np.random.default_rng(seed)
    observed = set(eligible_signature)
    observed_count = len(observed & target)
    null_rows: list[dict[str, object]] = []
    smds: list[float] = []
    for draw in range(1, draws + 1):
        available = np.ones(len(candidate_ids), dtype=bool)
        picked: list[str] = []
        for pool in neighbor_pools:
            feasible = pool[available[pool]][:10]
            if not len(feasible):
                feasible = np.flatnonzero(available)
            choice = int(rng.choice(feasible))
            available[choice] = False
            picked.append(str(candidate_ids[choice]))
        picked_set = set(picked)
        count = len(picked_set & target)
        control_mean = z.loc[picked].mean()
        smds.append(float((z.loc[eligible_signature].mean() - control_mean).abs().max()))
        null_rows.append({"layer": name, "draw": draw, "target_count": count, "target_rate": count / len(picked_set)})
    nulls = pd.DataFrame(null_rows)
    p = float((1 + nulls["target_count"].ge(observed_count).sum()) / (1 + len(nulls)))
    summary = {
        "analysis": "one-sided empirical matched-null enrichment; a higher target count is more enriched",
        "covariate_complete_signature_genes": len(eligible_signature),
        "covariate_complete_nonsignature_background_genes": len(candidate_ids),
        "observed_target_count": observed_count,
        "observed_target_rate": observed_count / len(eligible_signature),
        "null_median_target_count": float(nulls["target_count"].median()),
        "null_upper_95_target_count": float(nulls["target_count"].quantile(.95)),
        "matched_null_empirical_p": p,
        "max_absolute_standardized_mean_difference": {
            "median": float(np.median(smds)), "upper_95_percentile": float(np.quantile(smds, .95)), "maximum": float(np.max(smds)),
        },
    }
    return summary, nulls


def ids_from_first_column(path: Path) -> set[str]:
    header = pd.read_csv(path, sep="\t", compression="infer", nrows=0).columns
    column = "SEQ_ID" if "SEQ_ID" in header else header[0]
    table = pd.read_csv(path, sep="\t", compression="infer", usecols=[column], dtype=str)
    return {gene for gene in table[column].map(core) if gene}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--gff3", type=Path, required=True)
    parser.add_argument("--protein-evidence", type=Path, required=True)
    parser.add_argument("--rin-targets", type=Path, required=True)
    parser.add_argument("--occupancy-platform", type=Path, required=True)
    parser.add_argument("--differential", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--null-output", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260828)
    args = parser.parse_args()

    signature = set(pd.read_csv(args.signature, usecols=["gene_id"])["gene_id"].map(core)) - {""}
    covariates = matching_covariates(args.matrix, args.metadata, args.gff3)
    protein = pd.read_csv(args.protein_evidence)
    protein["gene_id"] = protein["gene_id"].map(core)
    protein_universe = set(protein.loc[truth(protein["protein_measured"]), "gene_id"]) - {""}
    protein_target = set(protein.loc[truth(protein["protein_differential_BR8_vs_MG"]), "gene_id"]) - {""}

    targets_raw = pd.read_csv(args.rin_targets)
    targets = targets_raw.copy()
    targets["gene_id"] = targets["gene_id"].map(core)
    admitted = set(targets.loc[truth(targets["admitted_binding"]), "gene_id"]) - {""}
    occupancy_universe = ids_from_first_column(args.occupancy_platform)
    # Keep the perturbation comparison in the same represented occupancy frame
    # as the frozen 2x2 table; two admitted source identifiers are absent from
    # the GPL15968 platform design and were not part of that original denominator.
    admitted &= occupancy_universe
    differential_raw = pd.read_csv(args.differential)
    # Preserve the exact identifier join used in the frozen GSE210589 summary
    # before canonicalising IDs for covariate matching.  This avoids admitting
    # version-suffixed DE rows that were not in its original tested denominator.
    admitted_raw = set(targets_raw.loc[truth(targets_raw["admitted_binding"]), "gene_id"])
    testable_raw = admitted_raw & set(differential_raw["gene_id"])
    differential = differential_raw.copy()
    differential["gene_id"] = differential["gene_id"].map(core)
    differential = differential.dropna(subset=["gene_id"])
    differential_target = {
        core(gene) for gene in differential_raw.loc[
            truth(differential_raw["differential"]) & differential_raw["gene_id"].isin(testable_raw), "gene_id"
        ] if core(gene)
    }
    perturbation_universe = admitted & {core(gene) for gene in testable_raw if core(gene)}

    layers = [
        ("protein_differential", protein_universe, protein_target),
        ("RIN_occupancy", occupancy_universe, admitted),
        ("rin1_differential_among_occupied_testable_genes", perturbation_universe, differential_target),
    ]
    output: dict[str, object] = {
        "scope": "Frozen-signature cross-layer sensitivity analyses matched for GSE267238 WT endpoint detection, absolute stage-effect magnitude, and Ensembl Plants release-63 gene length.",
        "frozen_signature_genes": len(signature), "draws": args.draws, "seed": args.seed,
        "layers": {},
    }
    all_nulls: list[pd.DataFrame] = []
    for index, (name, universe, target) in enumerate(layers):
        matched, nulls = match_layer(name, signature, universe, target, covariates, args.draws, args.seed + index)
        output["layers"][name] = {
            "assay_universe_genes": len(universe),
            "target_genes_in_assay_universe": len(target & universe),
            "unmatched_source_universe_2x2": raw_table(signature, universe, target),
            "matched_null": matched,
        }
        all_nulls.append(nulls)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.null_output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    pd.concat(all_nulls, ignore_index=True).to_csv(args.null_output, index=False, compression="infer")
    print(json.dumps(output, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
