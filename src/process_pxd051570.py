#!/usr/bin/env python3
"""Build gene-level proteome/phosphoproteome support from PXD051570 tables."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


SOLYC = re.compile(r"(Solyc\d+g\d+)")


def gene_id(values: pd.Series) -> pd.Series:
    return values.astype(str).str.extract(SOLYC, expand=False)


def fisher_greater(a: int, b: int, c: int, d: int) -> float:
    """Dependency-free one-sided Fisher exact p-value."""
    total, row1, col1 = a + b + c + d, a + b, a + c
    denominator = math.comb(total, row1)
    upper = min(row1, col1)
    return float(sum(math.comb(col1, x) * math.comb(total - col1, row1 - x) / denominator for x in range(a, upper + 1)))


def eligibility_matched_rate(matched_signature: set[str], background: set[str], target: set[str]) -> tuple[int, int, float]:
    """Calculate the observed target rate in the same eligibility frame as null draws."""
    measured = matched_signature & background
    if not measured:
        raise ValueError("No covariate-eligible signature genes are present in the measured background")
    target_count = len(measured & target)
    return len(measured), target_count, target_count / len(measured)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table-s1", type=Path, required=True)
    parser.add_argument("--table-s2", type=Path, required=True)
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--matched-nulls", type=Path, required=True)
    parser.add_argument("--gene-output", type=Path, required=True)
    parser.add_argument("--signature-output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260827)
    parser.add_argument("--permutations", type=int, default=1000)
    args = parser.parse_args()

    abundance = pd.read_excel(args.table_s1, sheet_name="Proteome normalized_full", header=1)
    abundance["gene_id"] = gene_id(abundance["Accession"])
    abundance = abundance.dropna(subset=["gene_id"])
    mg = [column for column in abundance if str(column).startswith("MG_")]
    ripe = [column for column in abundance if str(column).startswith("BR+8_")]
    if len(mg) != 3 or len(ripe) != 3:
        raise ValueError(f"Expected three MG and three BR+8 replicates; observed {mg} and {ripe}")
    abundance["protein_log2_BR8_vs_MG"] = np.log2(abundance[ripe].mean(axis=1) / abundance[mg].mean(axis=1))
    protein = abundance.groupby("gene_id", as_index=False).agg(
        protein_description=("Description", "first"),
        protein_log2_BR8_vs_MG=("protein_log2_BR8_vs_MG", "mean"),
    )
    protein["protein_measured"] = True

    de_protein = pd.read_excel(args.table_s2, sheet_name="DE_proteome", header=3)
    de_protein["gene_id"] = gene_id(de_protein["Accession"])
    de_protein = de_protein[(de_protein["Stage"] == "BR8/MG") & de_protein.gene_id.notna()].copy()
    de_protein["protein_de_log2_ratio"] = np.log2(pd.to_numeric(de_protein["Ratio"]))
    de_protein = de_protein.sort_values("Adjust P value").drop_duplicates("gene_id")
    de_protein = de_protein[["gene_id", "protein_de_log2_ratio", "P value", "Adjust P value"]].rename(
        columns={"P value": "protein_de_p", "Adjust P value": "protein_de_fdr"}
    )

    phospho = pd.read_excel(args.table_s2, sheet_name="DE_site-specific", header=1)
    phospho["gene_id"] = gene_id(phospho["Accession_protein"])
    phospho = phospho[(phospho["Group"] == "BR8/MG") & phospho.gene_id.notna()].copy()
    phospho["phospho_log2_ratio"] = np.log2(pd.to_numeric(phospho["Ratio"]))
    phospho["absolute_effect"] = phospho["phospho_log2_ratio"].abs()
    phospho = phospho.sort_values(["gene_id", "absolute_effect"], ascending=[True, False])
    phospho_gene = phospho.groupby("gene_id", as_index=False).agg(
        phosphosite_de_count=("Accession_phos", "size"),
        phosphosite_max_abs_log2_ratio=("absolute_effect", "max"),
        phosphosite_representative_log2_ratio=("phospho_log2_ratio", "first"),
        phosphosite_min_fdr=("Adjust P value", "min"),
        phosphosite_ids=("phos_site", lambda x: ";".join(sorted({str(v) for v in x if pd.notna(v) and str(v)}))),
    )

    evidence = protein.merge(de_protein, on="gene_id", how="left", validate="one_to_one").merge(
        phospho_gene, on="gene_id", how="outer", validate="one_to_one"
    )
    evidence["protein_measured"] = evidence["protein_measured"].fillna(False).astype(bool)
    evidence["protein_differential_BR8_vs_MG"] = evidence["protein_de_fdr"].notna()
    evidence["phosphosite_differential_BR8_vs_MG"] = evidence["phosphosite_de_count"].fillna(0).gt(0)

    signature = pd.read_csv(args.signature)
    signature_support = signature[["gene_id", "effect_source_a"]].merge(evidence, on="gene_id", how="left", validate="one_to_one")
    measured = signature_support[signature_support.protein_measured.eq(True)].copy()
    observed_direction = float((np.sign(measured.effect_source_a) == np.sign(measured.protein_log2_BR8_vs_MG)).mean())
    rng = np.random.default_rng(args.seed)
    randomized = np.array([
        (np.sign(measured.effect_source_a.to_numpy()) == np.sign(rng.permutation(measured.protein_log2_BR8_vs_MG.to_numpy()))).mean()
        for _ in range(args.permutations)
    ])

    background = set(evidence.loc[evidence.protein_measured, "gene_id"])
    protein_de = set(evidence.loc[evidence.protein_differential_BR8_vs_MG, "gene_id"]) & background
    phospho_de = set(evidence.loc[evidence.phosphosite_differential_BR8_vs_MG, "gene_id"]) & background
    sig = set(signature.gene_id) & background

    def overlap_stats(target: set[str]) -> dict[str, object]:
        a = len(sig & target)
        b = len(sig - target)
        c = len(target - sig)
        d = len(background - sig - target)
        odds = (a * d / (b * c)) if b and c else None
        return {"signature_overlap": a, "background_target_count": len(target), "odds_ratio": odds, "fisher_greater_p": fisher_greater(a, b, c, d)}

    nulls = pd.read_csv(args.matched_nulls, compression="infer")
    null_gene_column = next(
        (column for column in ("matched_null_gene_id", "matched_gene_id", "gene_id") if column in nulls),
        None,
    )
    if null_gene_column is None:
        raise ValueError("Matched-null table lacks a recognized gene identifier column")
    if "signature_gene_id" not in nulls:
        raise ValueError("Matched-null table lacks signature_gene_id, required for an eligibility-matched observed rate")
    matched_signature = set(nulls["signature_gene_id"].astype(str))
    if not matched_signature <= set(signature.gene_id.astype(str)):
        unexpected = sorted(matched_signature - set(signature.gene_id.astype(str)))[:5]
        raise ValueError(f"Matched-null signature genes are absent from the frozen signature: {unexpected}")
    def null_rates(target: set[str]) -> pd.Series:
        def rate(values: pd.Series) -> float:
            measured_genes = set(values) & background
            return len(measured_genes & target) / len(measured_genes) if measured_genes else np.nan
        return nulls.groupby("draw")[null_gene_column].apply(rate).dropna()

    protein_null = null_rates(protein_de)
    phospho_null = null_rates(phospho_de)
    protein_stats = overlap_stats(protein_de)
    phospho_stats = overlap_stats(phospho_de)
    matched_measured_count, matched_protein_count, protein_observed_rate = eligibility_matched_rate(matched_signature, background, protein_de)
    _, matched_phospho_count, phospho_observed_rate = eligibility_matched_rate(matched_signature, background, phospho_de)
    protein_stats["full_signature_measured_count"] = len(sig)
    phospho_stats["full_signature_measured_count"] = len(sig)
    protein_stats["matched_signature_eligible_count"] = len(matched_signature)
    phospho_stats["matched_signature_eligible_count"] = len(matched_signature)
    protein_stats["matched_signature_measured_count"] = matched_measured_count
    phospho_stats["matched_signature_measured_count"] = matched_measured_count
    protein_stats["matched_signature_target_count"] = matched_protein_count
    phospho_stats["matched_signature_target_count"] = matched_phospho_count
    protein_stats["signature_measured_rate"] = protein_observed_rate
    phospho_stats["signature_measured_rate"] = phospho_observed_rate
    protein_stats["matched_null_empirical_p"] = float((1 + (protein_null >= protein_observed_rate).sum()) / (1 + len(protein_null)))
    phospho_stats["matched_null_empirical_p"] = float((1 + (phospho_null >= phospho_observed_rate).sum()) / (1 + len(phospho_null)))

    summary = {
        "source": "PXD051570 / PMC13220761 Tables S1-S2",
        "context": "Ailsa Craig pericarp, BR+8 versus mature green, three biological replicates per stage",
        "quantified_protein_genes": len(background),
        "signature_protein_genes_measured": len(measured),
        "signature_protein_direction_concordance": observed_direction,
        "direction_permutation_p": float((1 + (randomized >= observed_direction).sum()) / (1 + len(randomized))),
        "protein_differential_enrichment": protein_stats,
        "phosphosite_differential_enrichment": phospho_stats,
        "independent_protein_support_gate": bool(
            (observed_direction > .5 and (1 + (randomized >= observed_direction).sum()) / (1 + len(randomized)) < .05)
            or protein_stats["matched_null_empirical_p"] < .05
            or phospho_stats["matched_null_empirical_p"] < .05
        ),
        "claim_boundary": "Protein abundance and phosphosite correspondence support cross-layer association, not transcriptional causality.",
    }
    args.gene_output.parent.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(args.gene_output, index=False)
    signature_support.to_csv(args.signature_output, index=False)
    args.summary.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
