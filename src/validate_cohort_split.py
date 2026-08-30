#!/usr/bin/env python3
"""Fail closed when a frozen external cohort enters development configuration."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def admitted(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"cohort_id", "independence_group", "analysis_role", "admission_status"}
    if not rows or required.difference(rows[0]):
        raise ValueError(f"{path} does not contain the required split columns")
    return [row for row in rows if row["admission_status"] == "eligible"]


def validate_split(development: list[dict[str, str]], external: list[dict[str, str]]) -> list[str]:
    errors: list[str] = []
    dev_ids = {row["cohort_id"] for row in development}
    ext_ids = {row["cohort_id"] for row in external}
    dev_groups = {row["independence_group"] for row in development}
    ext_groups = {row["independence_group"] for row in external}
    if dev_ids & ext_ids:
        errors.append(f"cohort ids occur in development and frozen external sets: {sorted(dev_ids & ext_ids)}")
    if dev_groups & ext_groups:
        errors.append(f"independence groups occur in development and frozen external sets: {sorted(dev_groups & ext_groups)}")
    frozen_roles = [row for row in development if "frozen_external" in row["analysis_role"]]
    if frozen_roles:
        errors.append(f"development config contains frozen external roles: {[row['cohort_id'] for row in frozen_roles]}")
    if not external or any("frozen_external" not in row["analysis_role"] for row in external):
        errors.append("external config must contain only explicitly frozen external roles")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--development", type=Path, required=True)
    parser.add_argument("--external", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    development = admitted(args.development)
    external = admitted(args.external)
    errors = validate_split(development, external)
    report = {
        "development_cohort_count": len(development),
        "development_independence_group_count": len({row["independence_group"] for row in development}),
        "frozen_external_cohort_count": len(external),
        "frozen_external_independence_group_count": len({row["independence_group"] for row in external}),
        "validation_errors": errors,
        "validation_passed": not errors,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
