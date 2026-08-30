#!/usr/bin/env python3
"""Generate gene-length, detection, and effect-magnitude matched null gene sets.

These are backgrounds for a later functional-enrichment test, not evidence of
enrichment by themselves.  Matching is performed within the independent
GSE267238 WT contrast so no annotation or pathway statistic enters selection.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


def core(value: str) -> str:
    hit = re.search(r"Solyc\d{2}g\d{6}", value)
    return hit.group(0) if hit else ""


def gene_lengths(path: Path) -> pd.Series:
    result: dict[str, int] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            ident = core(fields[8])
            if ident:
                result[ident] = abs(int(fields[4]) - int(fields[3])) + 1
    return pd.Series(result, name="gene_length_bp")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--gff3", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--n-sets", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260826)
    args = parser.parse_args()

    signature = pd.read_csv(args.signature)["gene_id"].map(core)
    matrix = pd.read_csv(args.matrix, compression="infer").set_index("gene_id")
    matrix.index = matrix.index.map(core)
    matrix = matrix[matrix.index != ""].groupby(level=0).mean()
    metadata = pd.read_csv(args.metadata)
    wt = metadata[metadata["genotype"] == "WT"]
    early = wt.loc[wt["stage_ordinal"] == 0, "sample_id"].tolist()
    late = wt.loc[wt["stage_ordinal"] == 1, "sample_id"].tolist()
    if len(early) < 2 or len(late) < 2:
        raise ValueError("GSE267238 WT endpoints need at least two biological libraries each.")
    cov = pd.DataFrame(index=matrix.index)
    cov["detection_rate"] = (matrix[early + late] > 1).mean(axis=1)
    cov["absolute_effect"] = (matrix[late].mean(axis=1) - matrix[early].mean(axis=1)).abs()
    cov = cov.join(gene_lengths(args.gff3), how="inner")
    cov = cov[cov["gene_length_bp"] > 0].copy()
    selected = sorted(set(signature).intersection(cov.index))
    if len(selected) < 400:
        raise ValueError(f"Only {len(selected)} signature genes have all matching covariates.")
    candidates = cov.drop(index=selected)
    values = cov[["detection_rate", "absolute_effect", "gene_length_bp"]].copy()
    values["gene_length_bp"] = np.log10(values["gene_length_bp"])
    mean, std = values.mean(), values.std().replace(0, 1)
    z = (values - mean) / std
    candidate_z = z.loc[candidates.index].to_numpy()
    selected_z = z.loc[selected].to_numpy()
    candidate_ids = candidates.index.to_numpy()
    neighbor_pools: list[np.ndarray] = []
    for vector in selected_z:
        distances = ((candidate_z - vector) ** 2).sum(axis=1)
        pool = np.argpartition(distances, 199)[:200]
        neighbor_pools.append(pool[np.argsort(distances[pool])])
    rng = np.random.default_rng(args.seed)
    rows: list[dict[str, object]] = []
    smds: list[float] = []
    for draw in range(args.n_sets):
        available = np.ones(len(candidate_ids), dtype=bool)
        picked: list[str] = []
        for gene, pool in zip(selected, neighbor_pools, strict=True):
            pool = pool[available[pool]][:10]
            if not len(pool):
                pool = np.flatnonzero(available)
            choice = int(rng.choice(pool))
            available[choice] = False
            picked.append(candidate_ids[choice])
            rows.append({"draw": draw + 1, "signature_gene_id": gene, "matched_null_gene_id": candidate_ids[choice]})
        control_mean = z.loc[picked].mean()
        smds.append(float((z.loc[selected].mean() - control_mean).abs().max()))
    output = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    report = {
        "scope": "Matched null gene sets for future functional-enrichment controls; no enrichment claim is calculated here.",
        "signature_genes_matched": len(selected),
        "background_genes": len(candidates),
        "null_sets": args.n_sets,
        "covariates": ["WT endpoint detection rate", "absolute WT breaker-minus-mature-green effect", "log10 Ensembl Plants release-63 gene length"],
        "matching_method": "random draw among the 10 closest available genes from a precomputed 200-neighbor pool in standardized three-covariate space",
        "max_absolute_standardized_mean_difference": {"median": float(np.median(smds)), "upper_95_percentile": float(np.percentile(smds, 95)), "maximum": float(np.max(smds))},
    }
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    main()
