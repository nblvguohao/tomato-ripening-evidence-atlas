#!/usr/bin/env python3
"""Validate the auditable public-study manifest before any analysis begins."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

REQUIRED_COLUMNS = {
    "study_id", "accession", "source_url", "organism", "assay", "fruit_or_tissue",
    "design", "intervention", "independence_group", "planned_role", "status", "eligibility_note",
}
ALLOWED_ROLES = {
    "train", "tune", "confirmation", "external_test", "frozen_external_test",
    "prospective_external", "sensitivity", "context",
    # Preserved negative: a study demoted by a preregistered audit gate rather
    # than dropped. GSE19326 recovered 156/464 frozen genes (33.6%) under strict
    # probe-sequence remapping, below the preregistered 70% coverage gate, and is
    # retained as a failed-candidate audit record.
    "failed_candidate",
}
ACCESSION_PATTERN = re.compile(r"^(?:(GSE|SRP)\d+|PMC\d+_MOESM\d+)$")


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("Manifest has no header row.")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    if not rows:
        raise ValueError("Manifest has no study rows.")
    return rows


def validate(rows: list[dict[str, str]]) -> list[str]:
    errors: list[str] = []
    for index, row in enumerate(rows, start=2):
        for column in REQUIRED_COLUMNS:
            if not row[column].strip():
                errors.append(f"row {index}: '{column}' is empty")
        if row["accession"] and not ACCESSION_PATTERN.match(row["accession"]):
            errors.append(f"row {index}: unsupported accession format '{row['accession']}'")
        if row["source_url"] and not row["source_url"].startswith("https://"):
            errors.append(f"row {index}: source_url must use HTTPS")
        if row["planned_role"] not in ALLOWED_ROLES:
            errors.append(f"row {index}: unknown planned_role '{row['planned_role']}'")

    for field in ("study_id", "accession"):
        duplicates = [value for value, count in Counter(row[field] for row in rows).items() if count > 1]
        if duplicates:
            errors.append(f"duplicate {field}: {', '.join(sorted(duplicates))}")

    grouped_roles: dict[str, set[str]] = {}
    for row in rows:
        grouped_roles.setdefault(row["independence_group"], set()).add(row["planned_role"])
    for group, roles in grouped_roles.items():
        if roles.intersection({"external_test", "frozen_external_test", "prospective_external"}) and roles.intersection({"train", "tune"}):
            errors.append(f"independence_group '{group}' spans external testing and development roles")
    return errors


def summary(rows: list[dict[str, str]]) -> dict[str, object]:
    return {
        "study_count": len(rows),
        "independent_group_count": len({row["independence_group"] for row in rows}),
        "roles": dict(Counter(row["planned_role"] for row in rows)),
        "status": dict(Counter(row["status"] for row in rows)),
        "accessions": [row["accession"] for row in rows],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    try:
        rows = load_manifest(args.manifest)
        errors = validate(rows)
    except (OSError, ValueError) as exc:
        print(f"MANIFEST INVALID: {exc}", file=sys.stderr)
        return 2

    report = summary(rows)
    report["validation_errors"] = errors
    rendered = json.dumps(report, ensure_ascii=False, indent=2)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    if errors:
        print("MANIFEST INVALID: correct the errors above before modelling.", file=sys.stderr)
        return 1
    print("MANIFEST VALID: source inventory and split roles are internally consistent.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
