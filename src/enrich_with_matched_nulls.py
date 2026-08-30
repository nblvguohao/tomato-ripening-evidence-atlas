#!/usr/bin/env python3
"""Empirical GO-term enrichment against precomputed matched null gene sets."""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


def core(value: str) -> str:
    return value.split(".")[0]


def bh(pvalues: pd.Series) -> pd.Series:
    order = pvalues.sort_values().index
    ranked = pvalues.loc[order].to_numpy() * len(pvalues) / np.arange(1, len(pvalues) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    output = pd.Series(index=order, data=np.minimum(ranked, 1.0))
    return output.reindex(pvalues.index)


def validate_matching_frame(full_signature: set[str], nulls: pd.DataFrame) -> tuple[set[str], pd.Series]:
    """Return the exact observed frame and validate one-to-one null draws."""
    required = {"draw", "signature_gene_id", "matched_null_gene_id"}
    missing = required - set(nulls.columns)
    if missing:
        raise ValueError(f"Matched-null table lacks required columns: {sorted(missing)}")
    matched_signature = set(nulls["signature_gene_id"])
    if not matched_signature <= full_signature:
        unexpected = sorted(matched_signature - full_signature)[:5]
        raise ValueError(f"Matched-null signature genes are absent from the frozen signature: {unexpected}")
    draw_sizes = nulls.groupby("draw")["matched_null_gene_id"].nunique()
    if draw_sizes.empty or not draw_sizes.eq(len(matched_signature)).all():
        raise ValueError(
            "Every matched-null draw must contain one unique control per eligible signature gene; "
            f"expected {len(matched_signature)}, observed {sorted(draw_sizes.unique())}"
        )
    return matched_signature, draw_sizes


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--annotations", type=Path, required=True)
    parser.add_argument("--matched-nulls", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-genes", type=int, default=5)
    args = parser.parse_args()
    full_signature = set(pd.read_csv(args.signature)["gene_id"].map(core))
    nulls = pd.read_csv(args.matched_nulls, compression="infer")
    nulls["signature_gene_id"] = nulls["signature_gene_id"].map(core)
    nulls["matched_null_gene_id"] = nulls["matched_null_gene_id"].map(core)
    matched_signature, draw_sizes = validate_matching_frame(full_signature, nulls)
    if args.annotations.name.endswith(".csv") or args.annotations.name.endswith(".csv.gz"):
        propagated = pd.read_csv(args.annotations, compression="infer", dtype=str).fillna("")
        if not {"gene_id", "propagated_terms"}.issubset(propagated.columns):
            raise ValueError("Propagated annotation CSV requires gene_id and propagated_terms columns")
        raw = propagated[["gene_id", "propagated_terms"]].rename(columns={"propagated_terms": "terms"})
        separator = ";"
    else:
        raw = pd.read_csv(args.annotations, sep="\t", header=None, names=["gene_id", "terms"], dtype=str).fillna("")
        separator = ","
    raw["gene_id"] = raw["gene_id"].map(core)
    ann = raw.assign(term=raw["terms"].str.split(separator)).explode("term")
    ann = ann[ann.term.str.startswith("GO:")][["gene_id", "term"]].drop_duplicates()
    universe = set(nulls["matched_null_gene_id"]) | matched_signature
    ann = ann[ann.gene_id.isin(universe)]
    # The observed set must use the same covariate-eligible signature frame as
    # every matched-null draw. The full frozen signature remains unchanged.
    observed = ann[ann.gene_id.isin(matched_signature)].groupby("term").gene_id.nunique()
    eligible = observed[observed >= args.min_genes]
    eligible_terms = set(eligible.index)
    annotation_map = ann.groupby("gene_id").term.apply(lambda values: set(values) & eligible_terms).to_dict()
    draw_ids = sorted(nulls["draw"].unique())
    null_counts = {term: np.zeros(len(draw_ids), dtype=int) for term in eligible_terms}
    for draw_index, (_, genes) in enumerate(nulls.groupby("draw", sort=True)):
        counts = Counter(term for gene in set(genes["matched_null_gene_id"]) for term in annotation_map.get(gene, set()))
        for term, count in counts.items():
            null_counts[term][draw_index] = count
    rows = []
    for term, count in eligible.items():
        draws = null_counts[term]
        p = (1 + int((draws >= count).sum())) / (1 + len(draws))
        rows.append({
            "go_id": term,
            "signature_gene_count": int(count),
            "matched_signature_gene_count": len(matched_signature),
            "full_frozen_signature_gene_count": len(full_signature),
            "null_draw_gene_count": int(draw_sizes.iloc[0]),
            "null_median_gene_count": float(np.median(draws)),
            "null_upper_95_gene_count": float(np.quantile(draws, .95)),
            "empirical_p": p,
        })
    result = pd.DataFrame(rows)
    if len(result):
        result["bh_fdr"] = bh(result["empirical_p"])
        result = result.sort_values(["bh_fdr", "empirical_p", "signature_gene_count"], ascending=[True, True, False])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"Tested {len(result)} GO terms with at least {args.min_genes} signature genes.")


if __name__ == "__main__":
    main()
