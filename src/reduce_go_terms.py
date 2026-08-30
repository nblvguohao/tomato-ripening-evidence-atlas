#!/usr/bin/env python3
"""Flag highly overlapping significant ancestor terms for transparent GO reduction."""

from __future__ import annotations

import argparse
from functools import lru_cache
from pathlib import Path

import pandas as pd

from propagate_go_annotations import parse_obo


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--enrichment", type=Path, required=True)
    parser.add_argument("--annotations", type=Path, required=True)
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--obo", type=Path, required=True)
    parser.add_argument("--term-table", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fdr", type=float, default=.05)
    parser.add_argument("--overlap", type=float, default=.8)
    args = parser.parse_args()

    enrichment = pd.read_csv(args.enrichment)
    parents, _ = parse_obo(args.obo)

    @lru_cache(maxsize=None)
    def ancestors(term: str) -> frozenset[str]:
        direct = parents.get(term, set())
        return frozenset(direct | {a for parent in direct for a in ancestors(parent)})

    signature = set(pd.read_csv(args.signature).gene_id)
    annotations = pd.read_csv(args.annotations, compression="infer", dtype=str).fillna("")
    long = annotations.assign(go_id=annotations.propagated_terms.str.split(";")).explode("go_id")
    long = long[long.gene_id.isin(signature) & long.go_id.str.startswith("GO:")]
    gene_sets = long.groupby("go_id").gene_id.apply(set).to_dict()
    significant = enrichment.loc[enrichment.bh_fdr < args.fdr, "go_id"].tolist()
    redundant_with: dict[str, str] = {}
    for ancestor in significant:
        ancestor_genes = gene_sets.get(ancestor, set())
        if not ancestor_genes:
            continue
        candidates = []
        for descendant in significant:
            if ancestor == descendant or ancestor not in ancestors(descendant):
                continue
            descendant_genes = gene_sets.get(descendant, set())
            overlap = len(ancestor_genes & descendant_genes) / len(ancestor_genes | descendant_genes) if descendant_genes else 0
            if overlap >= args.overlap:
                candidates.append((overlap, len(ancestors(descendant)), descendant))
        if candidates:
            redundant_with[ancestor] = max(candidates)[2]
    output = enrichment.copy()
    output["redundant_with_descendant"] = output.go_id.map(redundant_with).fillna("")
    output["representative_term"] = output["redundant_with_descendant"].eq("")
    names = pd.read_csv(args.term_table)
    output = output.merge(names, on="go_id", how="left", validate="one_to_one")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    print(f"Flagged {len(redundant_with)} of {len(significant)} significant propagated terms as highly overlapping ancestors.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
