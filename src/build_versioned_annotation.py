#!/usr/bin/env python3
"""Freeze a versioned annotation table for the frozen signature.

Stable Solyc locus cores are retained as analysis IDs.  A mapping row is kept
only when the supplied Ensembl export has exactly one matching gene record.
This avoids silently turning legacy probes or versioned transcript IDs into
ambiguous gene annotations.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


def locus_core(value: object) -> str:
    match = re.search(r"Solyc\d{2}g\d{6}", str(value))
    return match.group(0) if match else ""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--gene-metadata", type=Path, required=True)
    parser.add_argument("--ensembl-mapping", type=Path, required=True)
    parser.add_argument("--ensembl-release", default="63")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    signature = pd.read_csv(args.signature)
    metadata = pd.read_csv(args.gene_metadata, compression="infer")
    ensembl = pd.read_csv(args.ensembl_mapping, sep="\t", dtype=str).fillna("")
    required = {"Gene stable ID", "Gene name", "RefSeq mRNA ID"}
    missing = required.difference(ensembl.columns)
    if missing:
        raise ValueError(f"Ensembl mapping lacks columns: {sorted(missing)}")

    metadata = metadata.copy()
    metadata["gene_id"] = metadata["gene_id"].map(locus_core)
    metadata = metadata[metadata["gene_id"] != ""]
    metadata = metadata.drop_duplicates("gene_id", keep=False)

    ensembl = ensembl.copy()
    ensembl["gene_id"] = ensembl["Gene stable ID"].map(locus_core)
    ensembl = ensembl[ensembl["gene_id"] != ""]
    counts = ensembl.groupby("gene_id").size()
    ensembl = ensembl[ensembl["gene_id"].map(counts).eq(1)]
    ensembl = ensembl.rename(columns={
        "Gene stable ID": "ensembl_gene_stable_id",
        "Gene name": "ensembl_gene_name",
        "RefSeq mRNA ID": "refseq_mrna_id",
    })[["gene_id", "ensembl_gene_stable_id", "ensembl_gene_name", "refseq_mrna_id"]]

    annotation_columns = ["gene_id", "gene_name", "gene_biotype", "gene_description", "tf_family"]
    output = signature[["gene_id"]].merge(metadata[annotation_columns], on="gene_id", how="left", validate="one_to_one")
    output = output.merge(ensembl, on="gene_id", how="left", validate="one_to_one")
    output["annotation_source"] = "Ensembl Plants BioMart"
    output["annotation_release"] = str(args.ensembl_release)
    output["mapping_status"] = output["ensembl_gene_stable_id"].notna().map({True: "one_to_one", False: "not_mapped_in_export"})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    report = {
        "signature_genes": int(len(output)),
        "one_to_one_ensembl_mapped": int((output["mapping_status"] == "one_to_one").sum()),
        "not_mapped_in_export": int((output["mapping_status"] != "one_to_one").sum()),
        "ensembl_release": str(args.ensembl_release),
        "mapping_source_file": str(args.ensembl_mapping),
        "rule": "retain only one-to-one locus-core mappings from the supplied export",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
