#!/usr/bin/env python3
"""Prepare the frozen Micro-Tom GSE19326 validation matrix.

GPL4741's source annotation predates Solyc IDs.  This script retains only
probe sets supported by at least six exact 25-mer matches to one frozen-signature
canonical cDNA and no exact matches to another signature cDNA.  The resulting
mapping is intentionally scoped to the frozen signature and is not represented
as a whole-platform reannotation.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier = ""
    pieces: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier:
                    records[identifier] = "".join(pieces).upper()
                match = re.search(r"Solyc\d{2}g\d{6}", line)
                if not match:
                    raise ValueError(f"Cannot recover a Solyc locus from FASTA header: {line}")
                identifier = match.group(0)
                pieces = []
            else:
                pieces.append(line)
    if identifier:
        records[identifier] = "".join(pieces).upper()
    return records


def read_series(path: Path) -> tuple[list[str], pd.DataFrame]:
    titles: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = iter(handle)
        for line in lines:
            if line.startswith("!Sample_title"):
                titles = re.findall(r'"([^"]+)"', line)
            if line.startswith("!series_matrix_table_begin"):
                break
        header = next(lines).rstrip("\n").split("\t")
        rows: list[list[str]] = []
        for line in lines:
            if line.startswith("!series_matrix_table_end"):
                break
            rows.append(line.rstrip("\n").split("\t"))
    table = pd.DataFrame(rows, columns=[value.strip('"') for value in header])
    table = table.rename(columns={table.columns[0]: "probe_id"})
    table["probe_id"] = table["probe_id"].astype(str).str.strip('"')
    for column in table.columns[1:]:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    return titles, table


def frozen_microtom_metadata(sample_ids: list[str], titles: list[str]) -> pd.DataFrame:
    if len(sample_ids) != len(titles):
        raise ValueError("GSE19326 sample titles do not match matrix columns")
    records = []
    stage_order = {"mature_green": 1, "yellow": 2, "orange": 3, "red": 4}
    for sample_id, title in zip(sample_ids, titles, strict=True):
        value = title.lower()
        if "micro-tom flesh" not in value:
            continue
        stage = next((name for name, phrase in [
            ("mature_green", "mature green"), ("yellow", "yellow fruit"),
            ("orange", "orange fruit"), ("red", "red fruit"),
        ] if phrase in value), None)
        if stage is None:
            raise ValueError(f"Unable to determine GSE19326 stage: {title}")
        records.append({
            "sample_id": sample_id,
            "study_id": "GSE19326",
            "title": title,
            "cultivar": "Micro-Tom",
            "tissue": "flesh",
            "stage": stage,
            "stage_ordinal": stage_order[stage],
        })
    frame = pd.DataFrame(records)
    if frame.empty:
        raise ValueError("No Micro-Tom flesh samples were recovered")
    frame["replicate"] = frame.groupby("stage").cumcount() + 1
    counts = frame.groupby("stage").size().to_dict()
    if counts.get("mature_green", 0) < 2 or counts.get("red", 0) < 2:
        raise ValueError("GSE19326 requires at least two Micro-Tom flesh replicates at both endpoints")
    return frame


def map_signature_probes(probes: pd.DataFrame, transcripts: dict[str, str]) -> pd.DataFrame:
    required = {"probe_id", "sequence"}
    if missing := required.difference(probes.columns):
        raise ValueError(f"Probe table lacks columns: {sorted(missing)}")
    probes = probes[["probe_id", "sequence"]].dropna().copy()
    probes["sequence"] = probes["sequence"].astype(str).str.upper()
    probe_sequences = set(probes.sequence)
    probe_sequences.update(reverse_complement(sequence) for sequence in list(probe_sequences))
    hits: dict[str, set[str]] = defaultdict(set)
    for gene_id, transcript in transcripts.items():
        for index in range(max(0, len(transcript) - 24)):
            kmer = transcript[index : index + 25]
            if kmer in probe_sequences:
                hits[kmer].add(gene_id)
    probes["matched_genes"] = probes.sequence.map(lambda sequence: hits.get(sequence, set()).union(hits.get(reverse_complement(sequence), set())))
    records = []
    for probe_id, group in probes.groupby("probe_id", sort=False):
        by_gene: dict[str, int] = defaultdict(int)
        for genes in group.matched_genes:
            for gene_id in genes:
                by_gene[gene_id] += 1
        qualified = sorted(gene_id for gene_id, count in by_gene.items() if count >= 6)
        status = "one_to_one_signature_target" if len(qualified) == 1 else ("unmapped" if not qualified else "ambiguous_within_signature")
        records.append({
            "probe_id": probe_id,
            "gene_id": qualified[0] if status == "one_to_one_signature_target" else pd.NA,
            "mapping_status": status,
            "max_exact_probe_hits": max(by_gene.values(), default=0),
            "qualified_signature_targets": ";".join(qualified),
        })
    return pd.DataFrame(records)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--series", type=Path, required=True)
    parser.add_argument("--probes", type=Path, required=True)
    parser.add_argument("--signature-cdna", type=Path, required=True)
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()

    signature = pd.read_csv(args.signature)
    requested_genes = set(signature.gene_id)
    transcripts = {gene_id: sequence for gene_id, sequence in read_fasta(args.signature_cdna).items() if gene_id in requested_genes}
    if not transcripts:
        raise ValueError("No frozen-signature cDNA records were recovered")
    probes = pd.read_csv(args.probes, compression="infer")
    mapping = map_signature_probes(probes, transcripts)
    titles, expression = read_series(args.series)
    sample_ids = expression.columns[1:].tolist()
    metadata = frozen_microtom_metadata(sample_ids, titles)
    selected_ids = metadata.sample_id.tolist()
    joined = expression.merge(mapping, on="probe_id", how="left", validate="one_to_one")
    accepted = joined[joined.mapping_status == "one_to_one_signature_target"].dropna(subset=["gene_id"])
    gene_expression = accepted.groupby("gene_id", as_index=False)[selected_ids].median()
    rank01 = gene_expression[selected_ids].rank(axis=0, pct=True, method="average")
    output = pd.concat([gene_expression[["gene_id"]], rank01], axis=1)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(args.output_dir / "GSE19326_signature_probe_mapping.csv", index=False)
    output.to_csv(args.output_dir / "GSE19326_MicroTom_flesh_signature_rank01.csv.gz", index=False, compression="gzip")
    metadata.to_csv(args.output_dir / "GSE19326_MicroTom_flesh_sample_metadata.csv", index=False)
    stage_counts = metadata.groupby("stage").size().to_dict()
    audit = {
        "accession": "GSE19326",
        "platform": "GPL4741",
        "matrix_probe_rows": int(len(expression)),
        "matrix_sample_count": int(len(sample_ids)),
        "selected_microtom_flesh_samples": int(len(metadata)),
        "stage_sample_counts": {key: int(value) for key, value in stage_counts.items()},
        "signature_gene_count": int(len(signature)),
        "signature_cdna_available": int(len(transcripts)),
        "uniquely_mapped_signature_genes": int(len(output)),
        "signature_coverage": float(len(output) / len(signature)),
        "one_to_one_probe_sets": int((mapping.mapping_status == "one_to_one_signature_target").sum()),
        "mapping_rule": "at least 6 exact 25-mer probe hits to one frozen-signature canonical cDNA and no qualified hit to another signature cDNA",
        "mapping_scope": "unique within frozen-signature target cDNAs; not a complete GPL4741 reannotation",
        "admission_threshold_signature_coverage": 0.70,
        "admission": "pending_qc_and_preregistered_validation" if len(output) / len(signature) >= 0.70 else "ineligible_insufficient_signature_coverage",
    }
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    args.audit.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
