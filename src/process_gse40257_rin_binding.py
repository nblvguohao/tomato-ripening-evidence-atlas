#!/usr/bin/env python3
"""Create reproducible, source-bound RIN ChIP-chip binding evidence.

The GEO supplementary peak tables are analysed as supplied.  A gene is called
RIN-bound only when a peak with FDR <= 0.05 occurs in at least two independent
ChIP-chip replicates.  This program does *not* infer regulatory sign or
activation from occupancy.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


SOURCE_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40257"


def fisher_two_sided(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    """Dependency-free two-sided Fisher exact odds ratio and p-value."""
    total, row1, row2, col1 = a + b + c + d, a + b, c + d, a + c
    denominator = math.comb(total, row1)

    def probability(x: int) -> float:
        return math.comb(col1, x) * math.comb(total - col1, row1 - x) / denominator

    lower = max(0, row1 + col1 - total)
    upper = min(row1, col1)
    observed = probability(a)
    p_value = sum(probability(x) for x in range(lower, upper + 1) if probability(x) <= observed + 1e-15)
    odds_ratio = math.inf if b * c == 0 and a * d > 0 else (a * d / (b * c) if b * c else math.nan)
    return float(odds_ratio), float(min(p_value, 1.0))


def canonical_gene_id(values: pd.Series) -> pd.Series:
    """Return only stable ITAG-style Solyc identifiers from a source column."""
    return values.astype(str).str.extract(r"(Solyc\d{2}g\d{6})", expand=False)


def peak_targets(paths: list[Path], fdr_threshold: float = 0.05) -> pd.DataFrame:
    """Summarise replicate-level source-provided peaks by target gene."""
    records: list[pd.DataFrame] = []
    for path in paths:
        table = pd.read_csv(path, sep="\t", compression="infer")
        required = {"PEAK_FDR", "FEATURE_ATTR", "PEAK_SCORE", "PEAK_ID"}
        missing = required.difference(table.columns)
        if missing:
            raise ValueError(f"{path} lacks required peak columns: {sorted(missing)}")
        table = table.loc[pd.to_numeric(table["PEAK_FDR"], errors="coerce") <= fdr_threshold].copy()
        table["gene_id"] = canonical_gene_id(table["FEATURE_ATTR"])
        table = table.dropna(subset=["gene_id"])
        table["replicate"] = path.name.split("_")[0]
        records.append(table[["gene_id", "replicate", "PEAK_ID", "PEAK_SCORE", "PEAK_FDR"]])
    if not records:
        raise ValueError("At least one peak table is required")
    peaks = pd.concat(records, ignore_index=True)
    return peaks.groupby("gene_id", as_index=False).agg(
        binding_replicates=("replicate", "nunique"),
        source_peak_count=("PEAK_ID", "count"),
        max_peak_score=("PEAK_SCORE", "max"),
        min_peak_fdr=("PEAK_FDR", "min"),
        replicate_ids=("replicate", lambda items: ";".join(sorted(set(items)))),
    )


def array_gene_universe(path: Path) -> set[str]:
    """Return genes represented by the source tiling-array design."""
    design = pd.read_csv(path, sep="\t", compression="infer", usecols=["SEQ_ID"])
    return set(canonical_gene_id(design["SEQ_ID"]).dropna())


def enrichment(signature: set[str], bound: set[str], universe: set[str]) -> dict[str, float | int]:
    signature = signature.intersection(universe)
    bound = bound.intersection(universe)
    a = len(signature.intersection(bound))
    b = len(signature.difference(bound))
    c = len(bound.difference(signature))
    d = len(universe.difference(signature.union(bound)))
    odds_ratio, p_value = fisher_two_sided(a, b, c, d)
    return {
        "signature_genes_in_array_universe": len(signature),
        "bound_genes_in_array_universe": len(bound),
        "signature_bound_genes": a,
        "odds_ratio": float(odds_ratio),
        "fisher_p_value": float(p_value),
        "table_signature_bound": a,
        "table_signature_unbound": b,
        "table_nonsignature_bound": c,
        "table_nonsignature_unbound": d,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--peak-tables", type=Path, nargs="+", required=True)
    parser.add_argument("--array-design", type=Path, required=True)
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--consensus", type=Path, help="Four-cohort consensus used for graph-edge eligibility")
    parser.add_argument("--targets-output", type=Path, required=True)
    parser.add_argument("--evidence-output", type=Path, required=True)
    parser.add_argument("--binding-seeds-output", type=Path, help="Optional source-bound graph seed table")
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    targets = peak_targets(args.peak_tables)
    targets["admitted_binding"] = targets["binding_replicates"] >= 2
    targets = targets.sort_values(["admitted_binding", "binding_replicates", "max_peak_score"], ascending=[False, False, False])
    signature = set(pd.read_csv(args.signature, usecols=["gene_id"])["gene_id"])
    universe = array_gene_universe(args.array_design)
    admitted = set(targets.loc[targets["admitted_binding"], "gene_id"])
    target_signature = targets[targets["gene_id"].isin(signature)].copy()
    target_signature["evidence_layer"] = "chromatin"
    target_signature["assay"] = "RIN_ChIP_chip"
    target_signature["effect_direction"] = "binding_only"
    target_signature["context"] = "ripening_tomato_fruit;ITAG2;three_biological_replicates"
    target_signature["source_id"] = "GSE40257;PMID:23386264"
    target_signature["source_url"] = SOURCE_URL
    target_signature["biological_replicates"] = "three_source_reported"
    target_signature["admission_status"] = np.where(target_signature["admitted_binding"], "admitted", "below_two_replicate_threshold")
    target_signature["claim_boundary"] = "RIN occupancy only; no regulatory sign, causality, or activation magnitude inferred"
    target_signature["independent_perturbation_support"] = "false"
    evidence_columns = [
        "gene_id", "evidence_layer", "assay", "effect_direction", "context", "source_id", "source_url",
        "biological_replicates", "admission_status", "claim_boundary", "independent_perturbation_support",
        "binding_replicates", "source_peak_count", "max_peak_score", "min_peak_fdr", "replicate_ids",
    ]
    summary = {
        "source_id": "GSE40257",
        "source_url": SOURCE_URL,
        "assay": "RIN ChIP-chip",
        "source_biological_replicates": 3,
        "admission_rule": "PEAK_FDR <= 0.05 in at least 2 independent replicate peak-to-feature tables",
        "array_gene_universe": len(universe),
        "genes_with_any_source_peak": int(len(targets)),
        "genes_admitted_binding": int(targets["admitted_binding"].sum()),
        "signature_evidence_rows_including_below_threshold": int(len(target_signature)),
        "signature_genes_admitted_binding": int(target_signature["admitted_binding"].sum()),
        "signature_enrichment": enrichment(signature, admitted, universe),
        "claim_boundary": "Binding evidence is not assigned regulatory direction and cannot independently establish direct transcriptional regulation.",
    }
    for output in (args.targets_output, args.evidence_output, args.summary):
        output.parent.mkdir(parents=True, exist_ok=True)
    targets.to_csv(args.targets_output, index=False)
    target_signature[evidence_columns].to_csv(args.evidence_output, index=False)
    if args.consensus and args.binding_seeds_output:
        consensus = pd.read_csv(args.consensus, usecols=["gene_id", "all_four_same_direction"])
        graph_ids = set(consensus.loc[consensus["all_four_same_direction"].eq(True), "gene_id"])
        seeds = target_signature.loc[target_signature["admitted_binding"] & target_signature["gene_id"].isin(graph_ids)].copy()
        seeds = pd.DataFrame({
            "source_node": "RIN",
            "target_gene_id": seeds["gene_id"],
            "evidence_type": "in_vivo_promoter_binding",
            "direction": "binding_only",
            "context": "ripening_tomato_fruit;ITAG2;three_biological_replicates",
            "source_id": "GSE40257;PMID:23386264",
            "source_url": SOURCE_URL,
            "evidence_note": "Source peak-to-feature table: FDR <= 0.05 in at least two independent ChIP-chip replicates; occupancy only.",
            "admission_status": "admitted",
        })
        args.binding_seeds_output.parent.mkdir(parents=True, exist_ok=True)
        seeds.to_csv(args.binding_seeds_output, index=False)
        summary["four_cohort_consensus_genes_admitted_binding"] = int(len(seeds))
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
