#!/usr/bin/env python3
"""Validate provenance and claim boundaries for gene-level multi-omics evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


REQUIRED = {
    "gene_id", "evidence_layer", "assay", "effect_direction", "context", "source_id",
    "source_url", "biological_replicates", "admission_status", "claim_boundary",
    "independent_perturbation_support",
}
ALLOWED_LAYERS = {"protein", "phosphosite", "chromatin", "perturbation_expression"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    table = pd.read_csv(args.input, dtype=str).fillna("")
    errors = []
    missing = REQUIRED.difference(table.columns)
    if missing:
        errors.append(f"missing columns: {sorted(missing)}")
    else:
        for index, row in table.iterrows():
            row_number = index + 2
            for column in REQUIRED:
                if not row[column].strip():
                    errors.append(f"row {row_number}: empty {column}")
            if row["evidence_layer"] not in ALLOWED_LAYERS:
                errors.append(f"row {row_number}: unsupported evidence layer {row['evidence_layer']}")
            if not row["source_url"].startswith("https://"):
                errors.append(f"row {row_number}: source_url must use HTTPS")
            if row["independent_perturbation_support"].lower() not in {"true", "false"}:
                errors.append(f"row {row_number}: independent_perturbation_support must be true/false")
    report = {
        "record_count": len(table),
        "layers": table["evidence_layer"].value_counts().to_dict() if "evidence_layer" in table else {},
        "validation_errors": errors,
        "high_confidence_rule": "chromatin or protein evidence plus an independently sourced perturbation-expression record with compatible context",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
