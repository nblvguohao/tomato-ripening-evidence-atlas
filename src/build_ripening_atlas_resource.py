#!/usr/bin/env python3
"""Build a gene-level, evidence-tiered table for the ripening atlas resource."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--consensus", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--edges", type=Path, required=True)
    parser.add_argument("--studyshield", type=Path)
    parser.add_argument("--studyshield-significance", type=Path)
    parser.add_argument("--go-enrichment", type=Path)
    parser.add_argument("--go-annotations", type=Path)
    parser.add_argument("--protein-support", type=Path)
    parser.add_argument("--tissue-support", type=Path)
    parser.add_argument("--perturbation-support", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    consensus = pd.read_csv(args.consensus)
    annotation = pd.read_csv(args.annotation)
    edges = pd.read_csv(args.edges)
    regulatory = edges[edges["edge_type"] != "cross_cohort_expression_membership"].copy()
    regulatory = regulatory.rename(columns={"target": "gene_id", "source": "regulator"})
    regulatory = regulatory.groupby("gene_id", as_index=False).agg(
        published_regulator=("regulator", lambda values: ";".join(sorted(set(values)))),
        published_edge_type=("edge_type", lambda values: ";".join(sorted(set(values)))),
        published_edge_direction=("direction", lambda values: ";".join(sorted(set(values)))),
        published_edge_context=("context", lambda values: ";".join(sorted(set(values)))),
    )
    annotation = annotation.drop(columns=[column for column in ["gene_name", "gene_biotype", "gene_description", "tf_family"] if column in annotation])
    atlas = consensus.merge(annotation, on="gene_id", how="left", validate="one_to_one")
    atlas = atlas.merge(regulatory, on="gene_id", how="left", validate="one_to_one")
    if args.studyshield:
        shield = pd.read_csv(args.studyshield, compression="infer")
        atlas = atlas.merge(shield, on="gene_id", how="left", validate="one_to_one")
    if args.studyshield_significance:
        significance = pd.read_csv(args.studyshield_significance).drop(columns="studyshield_score", errors="ignore")
        significance = significance.rename(columns={"empirical_p": "studyshield_empirical_p", "bh_fdr": "studyshield_bh_fdr"})
        atlas = atlas.merge(significance, on="gene_id", how="left", validate="one_to_one")
    if args.go_enrichment and args.go_annotations:
        enriched = pd.read_csv(args.go_enrichment)
        if "representative_term" in enriched:
            enriched = enriched[enriched["representative_term"].eq(True)]
        enriched_ids = set(enriched.loc[enriched["bh_fdr"] < .05, "go_id"])
        if args.go_annotations.name.endswith(".csv") or args.go_annotations.name.endswith(".csv.gz"):
            raw_go = pd.read_csv(args.go_annotations, compression="infer", dtype=str).fillna("")
            raw_go = raw_go.rename(columns={"propagated_terms": "go_terms"})
            separator = ";"
        else:
            raw_go = pd.read_csv(args.go_annotations, sep="\t", header=None, names=["versioned_gene_id", "go_terms"], dtype=str).fillna("")
            raw_go["gene_id"] = raw_go["versioned_gene_id"].str.replace(r"\.\d+$", "", regex=True)
            separator = ","
        raw_go["significant_matched_null_go_terms"] = raw_go["go_terms"].map(
            lambda value: ";".join(sorted(enriched_ids.intersection(value.split(separator))))
        )
        go_by_gene = raw_go.groupby("gene_id", as_index=False)["significant_matched_null_go_terms"].agg(
            lambda values: ";".join(sorted({term for value in values for term in value.split(";") if term}))
        )
        atlas = atlas.merge(go_by_gene, on="gene_id", how="left", validate="one_to_one")
    atlas["atlas_claim"] = atlas["all_four_same_direction"].map({True: "four_cohort_expression_consensus", False: "frozen_signature_with_partial_replication"})
    atlas["regulatory_claim_boundary"] = atlas["published_edge_type"].notna().map({True: "source_bound_context_limited", False: "no_regulatory_claim"})
    if args.protein_support:
        protein = pd.read_csv(args.protein_support)
        protein_columns = [
            "gene_id", "protein_measured", "protein_log2_BR8_vs_MG", "protein_de_log2_ratio",
            "protein_de_p", "protein_de_fdr", "protein_differential_BR8_vs_MG",
            "phosphosite_de_count", "phosphosite_max_abs_log2_ratio",
            "phosphosite_representative_log2_ratio", "phosphosite_min_fdr",
            "phosphosite_ids", "phosphosite_differential_BR8_vs_MG",
        ]
        atlas = atlas.merge(protein[[column for column in protein_columns if column in protein]], on="gene_id", how="left", validate="one_to_one")
        atlas["protein_measured"] = atlas["protein_measured"].eq(True)
        protein_concordant = np.sign(atlas["effect_source_a"]) == np.sign(atlas["protein_log2_BR8_vs_MG"])
        atlas["protein_support_status"] = "not_measured"
        atlas.loc[atlas.protein_measured & protein_concordant, "protein_support_status"] = "measured_direction_concordant"
        atlas.loc[atlas.protein_measured & ~protein_concordant, "protein_support_status"] = "measured_direction_discordant"
        atlas.loc[atlas["protein_differential_BR8_vs_MG"].eq(True) & protein_concordant, "protein_support_status"] = "differential_direction_concordant"
        atlas["phosphosite_support_status"] = atlas["phosphosite_differential_BR8_vs_MG"].eq(True).map({True: "differential_site", False: "none"})
    else:
        atlas["protein_support_status"] = "not_assessed_public_matrix_pending"
        atlas["phosphosite_support_status"] = "not_assessed_public_matrix_pending"
    atlas["chromatin_support_status"] = "none"
    chromatin = atlas["published_edge_type"].fillna("").str.contains("chromatin|h3k27", case=False)
    binding = atlas["published_edge_type"].fillna("").str.contains("binding", case=False)
    atlas.loc[binding, "chromatin_support_status"] = "binding_only"
    atlas.loc[chromatin, "chromatin_support_status"] = "source_bound_chromatin_plus_expression_context"
    atlas["multiomics_high_confidence"] = False
    atlas["perturbation_support_status"] = "not_assessed"
    atlas["rin1_vs_wt_log2_fold_change"] = np.nan
    atlas["rin1_vs_wt_bh_fdr"] = np.nan
    atlas["regulatory_candidate_grade"] = "expression_only_or_binding_only"
    if args.perturbation_support:
        perturbation = pd.read_csv(args.perturbation_support)
        perturbation = perturbation[[
            "gene_id", "log2FoldChange", "padj", "differential", "evidence_grade",
        ]].rename(columns={
            "log2FoldChange": "rin1_vs_wt_log2_fold_change",
            "padj": "rin1_vs_wt_bh_fdr",
            "differential": "rin1_vs_wt_differential",
            "evidence_grade": "perturbation_support_status",
        })
        atlas = atlas.drop(columns=[
            "perturbation_support_status", "rin1_vs_wt_log2_fold_change", "rin1_vs_wt_bh_fdr",
        ])
        atlas = atlas.merge(perturbation, on="gene_id", how="left", validate="one_to_one")
        supported = atlas["perturbation_support_status"].eq("binding_plus_independent_perturbation_expression")
        atlas.loc[supported, "multiomics_high_confidence"] = True
        atlas.loc[supported, "regulatory_candidate_grade"] = "binding_plus_independent_perturbation_expression"
        atlas["perturbation_support_status"] = atlas["perturbation_support_status"].fillna("no_admitted_binding_record")
    if args.tissue_support:
        tissue = pd.read_csv(args.tissue_support)
        atlas = atlas.merge(tissue, on="gene_id", how="left", validate="one_to_one")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    atlas.to_csv(args.output, index=False)
    print(f"Wrote atlas resource with {len(atlas)} frozen signature genes and {len(regulatory)} source-bound regulatory records.")


if __name__ == "__main__":
    main()
