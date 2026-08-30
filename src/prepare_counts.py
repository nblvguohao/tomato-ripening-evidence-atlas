#!/usr/bin/env python3
"""Normalize source-specific public RNA-seq count matrices without pooling studies."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def clean_gene_id(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"(\.\d+)+$", "", regex=True)


def load_gse267238(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(path, sep="\t", compression="gzip")
    sample_columns = [column for column in raw.columns if column.count("_") == 2]
    metadata_columns = ["gene_id", "gene_name", "gene_biotype", "gene_description", "tf_family"]
    genes = raw[metadata_columns].copy()
    genes["gene_id"] = clean_gene_id(genes["gene_id"])
    counts = raw[["gene_id", *sample_columns]].copy()
    counts["gene_id"] = clean_gene_id(counts["gene_id"])
    counts = counts.groupby("gene_id", as_index=False)[sample_columns].sum()
    return counts, genes.drop_duplicates("gene_id")


def load_gse183836(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(path, compression="gzip")
    id_column = raw.columns[0]
    tomato = raw[raw[id_column].astype(str).str.startswith("Solyc")].copy()
    if tomato.empty:
        raise ValueError("No tomato Solyc identifiers found; host/pathogen filtering failed.")
    tomato = tomato.rename(columns={id_column: "gene_id"})
    tomato["gene_id"] = clean_gene_id(tomato["gene_id"])
    sample_columns = [column for column in tomato.columns if column.startswith("BB")]
    counts = tomato[["gene_id", *sample_columns]].groupby("gene_id", as_index=False).sum()
    genes = pd.DataFrame({"gene_id": counts["gene_id"], "source": "GSE183836_host_filtered"})
    return counts, genes


def load_gse235023(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(path, sep="\t", compression="gzip")
    raw = raw.rename(columns={raw.columns[0]: "gene_id"})
    raw["gene_id"] = clean_gene_id(raw["gene_id"])
    sample_columns = [column for column in raw.columns if column != "gene_id"]
    counts = raw[["gene_id", *sample_columns]].groupby("gene_id", as_index=False).sum()
    genes = pd.DataFrame({"gene_id": counts["gene_id"], "source": "GSE235023_public_counts"})
    return counts, genes


def load_gse285925(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Recover public integer counts while retaining only version-cleaned Solyc genes."""
    raw = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    tomato = raw[raw["gene_id"].astype(str).str.startswith("Solyc")].copy()
    tomato["gene_id"] = clean_gene_id(tomato["gene_id"])
    sample_columns = [column for column in tomato.columns if column.startswith("count.")]
    if len(sample_columns) != 12:
        raise ValueError(f"GSE285925: expected 12 count columns, found {len(sample_columns)}")
    counts = tomato[["gene_id", *sample_columns]].groupby("gene_id", as_index=False).sum()
    metadata_columns = [
        column for column in ("gene_id", "gene_name", "transcript_id", "GO", "Description", "trans_type")
        if column in tomato.columns
    ]
    genes = tomato[metadata_columns].drop_duplicates("gene_id")
    return counts, genes


def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    sample_columns = [column for column in counts.columns if column != "gene_id"]
    matrix = counts[sample_columns].astype(float)
    library_sizes = matrix.sum(axis=0)
    if (library_sizes <= 0).any():
        bad = ", ".join(library_sizes[library_sizes <= 0].index)
        raise ValueError(f"Non-positive library size for: {bad}")
    normalized = np.log2(matrix.div(library_sizes, axis=1) * 1_000_000 + 1)
    return pd.concat([counts[["gene_id"]], normalized], axis=1)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study", choices=("GSE267238", "GSE183836", "GSE235023", "GSE285925"), required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--sample-sheet", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    loaders = {
        "GSE267238": load_gse267238,
        "GSE183836": load_gse183836,
        "GSE235023": load_gse235023,
        "GSE285925": load_gse285925,
    }
    counts, gene_metadata = loaders[args.study](args.input)
    sample_sheet = pd.read_csv(args.sample_sheet)
    observed = set(counts.columns) - {"gene_id"}
    expected = set(sample_sheet["sample_id"])
    if not expected.issubset(observed):
        raise ValueError(f"Sample sheet contains absent matrix columns: {sorted(expected-observed)}")
    counts = counts[["gene_id", *sample_sheet["sample_id"].tolist()]]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    log_cpm(counts).to_csv(args.output_dir / f"{args.study}_log2cpm.csv.gz", index=False, compression="gzip")
    gene_metadata.to_csv(args.output_dir / f"{args.study}_gene_metadata.csv.gz", index=False, compression="gzip")
    sample_sheet.to_csv(args.output_dir / f"{args.study}_sample_metadata.csv", index=False)
    print(f"Prepared {args.study}: {len(counts):,} genes and {len(expected)} selected samples.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
