#!/usr/bin/env python3
"""Attach versioned public gene descriptions to a frozen signature table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--gene-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    signature = pd.read_csv(args.signature)
    metadata = pd.read_csv(args.gene_metadata, compression="gzip")
    annotation_columns = [column for column in ["gene_id", "gene_name", "gene_biotype", "gene_description", "tf_family"] if column in metadata]
    output = signature.merge(metadata[annotation_columns].drop_duplicates("gene_id"), on="gene_id", how="left", validate="one_to_one")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    mapped = output["gene_description"].notna().sum() if "gene_description" in output else 0
    print(f"Annotated {mapped}/{len(output)} signature genes using the supplied source metadata.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
