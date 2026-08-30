#!/usr/bin/env python3
"""Validate all frozen GSE210589 inputs and per-sample processing gates."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--quant-root", type=Path, required=True)
    parser.add_argument("--fastp-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--minimum-retention", type=float, default=0.90)
    parser.add_argument("--minimum-mapping", type=float, default=70.0)
    args = parser.parse_args()

    samples = pd.read_csv(args.samples, dtype=str)
    required = {"sample_id", "sra_run", "condition", "biological_replicate", "inclusion_status"}
    missing = required.difference(samples.columns)
    if missing:
        raise ValueError(f"sample sheet is missing required columns: {sorted(missing)}")
    samples = samples[samples["inclusion_status"].eq("admitted")].copy()
    errors: list[str] = []
    if len(samples) != 6 or samples["sra_run"].nunique() != 6:
        errors.append("expected six unique admitted runs")
    group_sizes = samples.groupby("condition")["biological_replicate"].nunique().to_dict()
    if group_sizes != {"WT_34DPA": 3, "rin1_34DPA": 3}:
        errors.append(f"unexpected biological replicate counts: {group_sizes}")

    records = []
    for row in samples.itertuples(index=False):
        run_id = row.sra_run
        fastp_path = args.fastp_root / run_id / "fastp.json"
        meta_path = args.quant_root / run_id / "aux_info" / "meta_info.json"
        complete_path = args.quant_root / run_id / "COMPLETE"
        if not fastp_path.exists() or not meta_path.exists() or not complete_path.exists():
            errors.append(f"{run_id}: missing fastp, Salmon metadata, or COMPLETE marker")
            continue
        fastp = json.loads(fastp_path.read_text())
        meta = json.loads(meta_path.read_text())
        before = int(fastp["summary"]["before_filtering"]["total_reads"])
        after = int(fastp["summary"]["after_filtering"]["total_reads"])
        retention = after / before if before else 0.0
        mapping = float(meta["percent_mapped"])
        if retention < args.minimum_retention:
            errors.append(f"{run_id}: read retention {retention:.4f} below {args.minimum_retention:.4f}")
        if mapping < args.minimum_mapping:
            errors.append(f"{run_id}: mapping {mapping:.2f} below {args.minimum_mapping:.2f}")
        records.append({
            "sample_id": row.sample_id,
            "sra_run": run_id,
            "condition": row.condition,
            "biological_replicate": int(row.biological_replicate),
            "reads_before_filtering": before,
            "reads_after_filtering": after,
            "read_retention": retention,
            "salmon_percent_mapped": mapping,
            "complete": True,
        })
    report = {
        "source_id": "GSE210589",
        "expected_runs": 6,
        "validated_runs": len(records),
        "biological_replicates_by_condition": group_sizes,
        "minimum_read_retention": min((row["read_retention"] for row in records), default=None),
        "minimum_salmon_percent_mapped": min((row["salmon_percent_mapped"] for row in records), default=None),
        "sample_qc": records,
        "validation_errors": errors,
        "validation_passed": not errors,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
