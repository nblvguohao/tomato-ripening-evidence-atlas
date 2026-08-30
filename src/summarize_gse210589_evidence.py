#!/usr/bin/env python3
"""Add reproducible QC and admitted-target enrichment to GSE210589 summary."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from process_pxd051570 import fisher_greater


def enrichment_table(
    differential: pd.DataFrame,
    admitted_targets: pd.DataFrame,
    signature_binding: pd.DataFrame,
) -> tuple[list[list[int]], int]:
    admitted = admitted_targets[admitted_targets["admitted_binding"].eq(True)][["gene_id"]].drop_duplicates()
    tested = admitted.merge(differential[["gene_id", "differential"]], on="gene_id", how="inner")
    signature_ids = set(signature_binding["gene_id"])
    is_signature = tested["gene_id"].isin(signature_ids)
    signature_de = int(tested.loc[is_signature, "differential"].eq(True).sum())
    signature_not_de = int(is_signature.sum() - signature_de)
    other_de = int(tested.loc[~is_signature, "differential"].eq(True).sum())
    other_not_de = int((~is_signature).sum() - other_de)
    return [[signature_de, signature_not_de], [other_de, other_not_de]], len(tested)


def pair_values(matrix: pd.DataFrame, first: list[str], second: list[str], within: bool) -> list[float]:
    values: list[float] = []
    for i, left in enumerate(first):
        for j, right in enumerate(second):
            if within and j <= i:
                continue
            values.append(float(matrix.loc[left, right]))
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--processing-audit", type=Path, required=True)
    parser.add_argument("--differential", type=Path, required=True)
    parser.add_argument("--rin-targets", type=Path, required=True)
    parser.add_argument("--signature-binding", type=Path, required=True)
    parser.add_argument("--sample-correlations", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    args = parser.parse_args()

    summary = json.loads(args.summary.read_text())
    audit = json.loads(args.processing_audit.read_text())
    differential = pd.read_csv(args.differential)
    targets = pd.read_csv(args.rin_targets)
    signature = pd.read_csv(args.signature_binding)
    table, target_count = enrichment_table(differential, targets, signature)
    a, b = table[0]
    c, d = table[1]
    odds_ratio = float("inf") if b * c == 0 and a * d > 0 else (a * d / (b * c) if b * c else float("nan"))
    p_value = fisher_greater(a, b, c, d)

    correlations = pd.read_csv(args.sample_correlations).set_index("sample_id")
    samples = pd.read_csv(args.samples)
    wt = samples.loc[samples.condition.eq("WT_34DPA"), "sra_run"].tolist()
    rin = samples.loc[samples.condition.eq("rin1_34DPA"), "sra_run"].tolist()
    within_wt = pair_values(correlations, wt, wt, within=True)
    within_rin = pair_values(correlations, rin, rin, within=True)
    between = pair_values(correlations, wt, rin, within=False)

    high_grade = signature[signature["evidence_grade"].eq("binding_plus_independent_perturbation_expression")]
    summary.update({
        "processing_validation_passed": bool(audit["validation_passed"]),
        "minimum_read_retention": float(audit["minimum_read_retention"]),
        "minimum_salmon_percent_mapped": float(audit["minimum_salmon_percent_mapped"]),
        "within_WT_spearman_range": [min(within_wt), max(within_wt)],
        "within_rin1_spearman_range": [min(within_rin), max(within_rin)],
        "between_group_spearman_range": [min(between), max(between)],
        "admitted_RIN_targets_tested": int(target_count),
        "signature_RIN_bound_differential": int(table[0][0]),
        "signature_RIN_bound_not_differential": int(table[0][1]),
        "other_RIN_bound_differential": int(table[1][0]),
        "other_RIN_bound_not_differential": int(table[1][1]),
        "signature_vs_other_RIN_target_differential_odds_ratio": float(odds_ratio),
        "signature_vs_other_RIN_target_differential_fisher_p": float(p_value),
        "perturbation_statistical_support_gate": bool(odds_ratio > 1 and p_value < 0.05),
        "high_grade_lower_in_rin1": int(high_grade["log2FoldChange"].lt(0).sum()),
        "high_grade_higher_in_rin1": int(high_grade["log2FoldChange"].gt(0).sum()),
    })
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
