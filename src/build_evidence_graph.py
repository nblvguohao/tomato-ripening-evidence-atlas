#!/usr/bin/env python3
"""Export a transparent evidence graph from consensus genes and source-bound edges.

Expression-consensus edges are membership edges, not regulatory edges. Published
binding edges retain their source and do not acquire a direction unless the
primary source explicitly establishes one.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--consensus", type=Path, required=True)
    parser.add_argument("--binding-seeds", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    consensus = pd.read_csv(args.consensus)
    replicated = consensus[consensus["all_four_same_direction"]].copy()
    signature_node = pd.DataFrame([{
        "node_id": "four_cohort_expression_signature",
        "label": "Four-cohort expression consensus",
        "node_type": "evidence_set",
        "candidate_rank": 0,
        "gene_description": "",
        "tf_family": "",
    }])
    genes = replicated[["gene_id", "gene_description", "tf_family", "candidate_rank"]].rename(columns={"gene_id": "node_id"})
    genes["label"] = genes["node_id"]
    genes["node_type"] = "gene"
    genes = genes[["node_id", "label", "node_type", "candidate_rank", "gene_description", "tf_family"]]

    seeds = pd.read_csv(args.binding_seeds)
    required_seed_columns = {
        "source_node", "target_gene_id", "evidence_type", "direction", "context",
        "source_id", "source_url", "evidence_note", "admission_status",
    }
    missing_seed_columns = required_seed_columns.difference(seeds.columns)
    if missing_seed_columns:
        raise ValueError(f"Binding seed table is missing columns: {sorted(missing_seed_columns)}")
    unknown_targets = set(seeds["target_gene_id"]).difference(set(genes["node_id"]))
    if unknown_targets:
        raise ValueError(f"Binding seed targets are not four-cohort consensus genes: {sorted(unknown_targets)}")
    regulators = pd.DataFrame({"node_id": sorted(set(seeds["source_node"]))})
    regulators["label"] = regulators["node_id"]
    regulators["node_type"] = "regulator_label"
    regulators["candidate_rank"] = 0
    regulators["gene_description"] = ""
    regulators["tf_family"] = ""
    nodes = pd.concat([signature_node, genes, regulators], ignore_index=True).drop_duplicates("node_id")

    membership_edges = pd.DataFrame({
        "source": "four_cohort_expression_signature",
        "target": replicated["gene_id"],
        "edge_type": "cross_cohort_expression_membership",
        "direction": "not_regulatory",
        "context": "four_expression_sources_same_direction",
        "source_id": "this_project",
        "source_url": pd.NA,
        "evidence_note": "Membership requires same direction in source A, confirmation source B, GSE267238, and GSE128739.",
        "admission_status": "admitted",
    })
    binding_edges = seeds.rename(columns={
        "source_node": "source",
        "target_gene_id": "target",
        "evidence_type": "edge_type",
    })
    edges = pd.concat([membership_edges, binding_edges], ignore_index=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    nodes.to_csv(args.output_dir / "evidence_graph_nodes.csv", index=False)
    edges.to_csv(args.output_dir / "evidence_graph_edges.csv", index=False)
    print(f"Wrote {len(nodes)} nodes and {len(edges)} edges; {len(seeds)} edge(s) are source-bound published binding evidence.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
