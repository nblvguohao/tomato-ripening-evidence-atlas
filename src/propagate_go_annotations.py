#!/usr/bin/env python3
"""Propagate direct ITAG GO terms through a versioned OBO is_a/part_of hierarchy."""

from __future__ import annotations

import argparse
from functools import lru_cache
from pathlib import Path

import pandas as pd


def parse_obo(path: Path) -> tuple[dict[str, set[str]], dict[str, tuple[str, str]]]:
    parents: dict[str, set[str]] = {}
    terms: dict[str, tuple[str, str]] = {}
    current: dict[str, object] | None = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line == "[Term]":
            if current and current.get("id"):
                identifier = str(current["id"])
                parents[identifier] = set(current.get("parents", set()))
                terms[identifier] = (str(current.get("name", "")), str(current.get("namespace", "")))
            current = {"parents": set()}
        elif current is not None and line.startswith("id: GO:"):
            current["id"] = line.split("id: ", 1)[1]
        elif current is not None and line.startswith("name: "):
            current["name"] = line.split("name: ", 1)[1]
        elif current is not None and line.startswith("namespace: "):
            current["namespace"] = line.split("namespace: ", 1)[1]
        elif current is not None and line.startswith("is_a: GO:"):
            current["parents"].add(line.split()[1])
        elif current is not None and line.startswith("relationship: part_of GO:"):
            current["parents"].add(line.split()[2])
    if current and current.get("id"):
        identifier = str(current["id"])
        parents[identifier] = set(current.get("parents", set()))
        terms[identifier] = (str(current.get("name", "")), str(current.get("namespace", "")))
    return parents, terms


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotations", type=Path, required=True)
    parser.add_argument("--obo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--term-table", type=Path, required=True)
    args = parser.parse_args()
    if not args.obo.is_file() or args.obo.stat().st_size == 0:
        raise FileNotFoundError("A non-empty, versioned GO OBO file is required")
    parents, terms = parse_obo(args.obo)

    @lru_cache(maxsize=None)
    def ancestors(term: str) -> frozenset[str]:
        direct = parents.get(term, set())
        return frozenset(direct | {ancestor for parent in direct for ancestor in ancestors(parent)})

    raw = pd.read_csv(args.annotations, sep="\t", header=None, names=["versioned_gene_id", "direct_terms"], dtype=str).fillna("")
    raw["gene_id"] = raw["versioned_gene_id"].str.replace(r"\.\d+$", "", regex=True)
    raw["propagated_terms"] = raw["direct_terms"].map(lambda value: ";".join(sorted({expanded for term in value.split(",") if term for expanded in ({term} | set(ancestors(term)))})))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    raw[["gene_id", "versioned_gene_id", "direct_terms", "propagated_terms"]].to_csv(args.output, index=False, compression="gzip" if args.output.suffix == ".gz" else None)
    term_table = pd.DataFrame([{"go_id": identifier, "go_name": value[0], "go_namespace": value[1]} for identifier, value in terms.items()])
    term_table.to_csv(args.term_table, index=False)
    print(f"Propagated annotations for {len(raw)} genes through {len(terms)} GO terms.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
