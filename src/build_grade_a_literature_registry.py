#!/usr/bin/env python3
"""Create a complete, status-explicit literature-screen registry for Grade A genes."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--existing", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    candidates = pd.read_csv(args.candidates)
    grade_a = candidates[candidates.candidate_grade.eq("A")].copy()
    if len(grade_a) != 18:
        raise ValueError(f"Expected 18 Grade A candidates, found {len(grade_a)}")
    existing = pd.read_csv(args.existing, dtype=str).fillna("")
    existing = existing.drop_duplicates("gene_id", keep="last")
    columns = [
        "gene_id", "display_name", "candidate_grade", "submission_priority_rank",
        "literature_screen_status", "literature_screen_scope", "literature_screen_source",
        "literature_screen_url", "literature_screen_conclusion", "screen_date",
    ]
    output = grade_a[["gene_id", "display_name", "candidate_grade", "submission_priority_rank"]].merge(
        existing.drop(columns=[column for column in ["display_name", "candidate_grade", "submission_priority_rank"] if column in existing], errors="ignore"),
        on="gene_id", how="left", validate="one_to_one",
    )
    for column in columns[4:]:
        if column not in output:
            output[column] = ""
        output[column] = output[column].fillna("")
    pending = output.literature_screen_status.eq("")
    output.loc[pending, "literature_screen_status"] = "pending_full_alias_and_full_text_screen"
    output.loc[pending, "literature_screen_scope"] = (
        "Pending: exact Solyc ID, aliases, gene symbol, tomato fruit-ripening context, and primary full-text review"
    )
    output.loc[pending, "literature_screen_conclusion"] = (
        "No manuscript novelty or direct-regulation claim is permitted until this row is completed."
    )
    output = output[columns].sort_values("submission_priority_rank")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    print(f"Wrote {len(output)} Grade A literature-screen rows; pending={int(pending.sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
