#!/usr/bin/env python3
"""Export a fixed, cross-study concordant mature-green-to-ripe signature."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from compare_stage_effects import load


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-a-matrix", type=Path, required=True)
    parser.add_argument("--study-a-metadata", type=Path, required=True)
    parser.add_argument("--study-b-matrix", type=Path, required=True)
    parser.add_argument("--study-b-metadata", type=Path, required=True)
    parser.add_argument("--top-k", type=int, default=500)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    effect_a = load(args.study_a_matrix, args.study_a_metadata, early=0, late=3)
    effect_b = load(args.study_b_matrix, args.study_b_metadata, early=0, late=3)
    common = effect_a.index.intersection(effect_b.index).sort_values()
    table = pd.DataFrame({"gene_id": common, "effect_source_a": effect_a.loc[common].to_numpy(), "effect_source_b": effect_b.loc[common].to_numpy()})
    table["same_direction"] = (table["effect_source_a"] * table["effect_source_b"]) > 0
    table["minimum_absolute_effect"] = table[["effect_source_a", "effect_source_b"]].abs().min(axis=1)
    selected = table.loc[table["effect_source_a"].abs().nlargest(args.top_k).index].copy()
    selected = selected[selected["same_direction"]].sort_values("minimum_absolute_effect", ascending=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output, index=False)
    print(f"Wrote {len(selected)} same-direction genes from a locked top-{args.top_k} source-A signature.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
