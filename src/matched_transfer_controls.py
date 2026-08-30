#!/usr/bin/env python3
"""Evaluate frozen-signature transfer against effect-matched control gene sets.

The frozen signature is the 464 genes with the largest absolute GSE42783
red-ripe-minus-mature-green effect that also retained direction in GSE108415.
Every admitted transfer cohort is itself a tomato fruit ripening contrast, so
high direction agreement cannot be interpreted without a control that shares
the selection variable.  The permutation null used in the frozen validation
files shuffles imported effects *within* the measured signature genes; it tests
pairing, not transferability relative to comparable genes.

This script adds the two controls the frozen outputs lack:

  * matched draws - 1,000 non-signature gene sets of the same size, matched on
    the selection variable (absolute source effect), source expression level,
    and log10 gene length;
  * rank 501-1000 - the next 500 genes by absolute source effect, a fixed,
    parameter-free control that asks whether transferability is a property of
    this set or a monotone function of source effect magnitude.

For each cohort and each control the script recomputes exactly the statistics
the manuscript reports for the signature: direction agreement against the sign
of the source effect, and the Spearman correlation of effect ranks.

Nothing here changes signature membership, reference direction, cohort roles,
or any admission gate.  It is a read-only sensitivity analysis over frozen
inputs.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

SOLYC = re.compile(r"(Solyc\d{2}g\d{6})")


def core(value: object) -> str:
    match = SOLYC.search(str(value))
    return match.group(1) if match else ""


def gene_lengths(path: Path) -> pd.Series:
    lengths: dict[str, int] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            gene = core(fields[8])
            if gene:
                lengths[gene] = abs(int(fields[4]) - int(fields[3])) + 1
    return pd.Series(lengths, name="gene_length_bp")


def load_matrix(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, compression="infer").set_index("gene_id")
    matrix.index = matrix.index.map(core)
    return matrix[matrix.index != ""].groupby(level=0).mean()


def endpoint_effect(matrix: pd.DataFrame, metadata: pd.DataFrame, genotype: str,
                    early: int, late: int) -> pd.Series:
    rows = metadata[metadata["genotype"] == genotype]
    early_ids = rows.loc[rows["stage_ordinal"] == early, "sample_id"].tolist()
    late_ids = rows.loc[rows["stage_ordinal"] == late, "sample_id"].tolist()
    early_ids = [s for s in early_ids if s in matrix.columns]
    late_ids = [s for s in late_ids if s in matrix.columns]
    if len(early_ids) < 2 or len(late_ids) < 2:
        raise ValueError("each endpoint needs at least two biological libraries")
    return matrix[late_ids].mean(axis=1) - matrix[early_ids].mean(axis=1)


def transfer_statistics(source: pd.Series, imported: pd.Series, genes: list[str]) -> dict[str, float] | None:
    """Direction agreement and effect-rank correlation, as the manuscript defines them."""
    usable = [g for g in genes if g in source.index and g in imported.index]
    if len(usable) < 30:
        return None
    a = source.loc[usable]
    b = imported.loc[usable]
    return {
        "genes": len(usable),
        "direction": float((np.sign(a) == np.sign(b)).mean()),
        "rho": float(a.rank().corr(b.rank(), method="pearson")),
    }


def matched_pools(selected: list[str], candidates: np.ndarray, z: pd.DataFrame,
                  pool_size: int = 200) -> list[np.ndarray]:
    candidate_z = z.loc[candidates].to_numpy()
    pools = []
    size = min(pool_size, len(candidates))
    for vector in z.loc[selected].to_numpy():
        distances = ((candidate_z - vector) ** 2).sum(axis=1)
        pool = np.argpartition(distances, size - 1)[:size]
        pools.append(pool[np.argsort(distances[pool])])
    return pools


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--source-matrix", type=Path, required=True)
    parser.add_argument("--source-metadata", type=Path, required=True)
    parser.add_argument("--gff3", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--null-output", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260828)
    args = parser.parse_args()

    signature = [core(g) for g in pd.read_csv(args.signature)["gene_id"]]
    signature = sorted({g for g in signature if g})

    # Source effect, recomputed for every gene rather than only for the frozen set.
    source_matrix = load_matrix(args.source_matrix)
    source_metadata = pd.read_csv(args.source_metadata)
    mature_green = source_metadata.loc[source_metadata["stage"] == "mature_green", "sample_id"].tolist()
    red_ripe = source_metadata.loc[source_metadata["stage"] == "red_ripe", "sample_id"].tolist()
    source_effect = source_matrix[red_ripe].mean(axis=1) - source_matrix[mature_green].mean(axis=1)
    source_level = source_matrix[mature_green + red_ripe].mean(axis=1)

    frozen = pd.read_csv(args.signature).set_index("gene_id")["effect_source_a"]
    shared = source_effect.index.intersection(frozen.index)
    drift = float((source_effect.loc[shared] - frozen.loc[shared]).abs().max())
    if drift > 1e-9:
        raise ValueError(f"recomputed source effect disagrees with the frozen file (max {drift})")

    covariates = pd.DataFrame({
        "absolute_source_effect": source_effect.abs(),
        "source_level": source_level,
    }).join(gene_lengths(args.gff3), how="inner")
    covariates = covariates[covariates["gene_length_bp"] > 0]
    covariates["gene_length_bp"] = np.log10(covariates["gene_length_bp"])
    z = (covariates - covariates.mean()) / covariates.std().replace(0, 1)

    # rank 501-1000 by absolute source effect: a fixed, parameter-free control
    ordered = source_effect.abs().sort_values(ascending=False)
    rank_control = [g for g in ordered.index[500:1000] if g not in set(signature)]

    cohorts = [
        ("GSE267238", "results/signature_validation_GSE267238.json",
         "data/processed/GSE267238_log2cpm.csv.gz", "data/processed/GSE267238_sample_metadata.csv",
         "WT", 0, 1, "admitted"),
        # day 53 minus day 41, matching the frozen contrast definition
        ("GSE128739", "results/signature_validation_GSE128739.json",
         "data/processed/GSE128739_log2cpm.csv.gz", "data/processed/GSE128739_sample_metadata.csv",
         "Moneymaker", 0, 3, "admitted"),
        ("GSE235023", "results/signature_validation_GSE235023.json",
         "data/processed/GSE235023_log2cpm.csv.gz", "data/processed/GSE235023_sample_metadata.csv",
         "AC_exocarp", 0, 1, "admitted"),
        ("GSE78733", None,
         "data/processed/GSE78733_rank01.csv.gz", "data/processed/GSE78733_sample_metadata.csv",
         "WT", 1, 3, "untouched"),
        ("GSE285925_No19", "results/signature_validation_GSE285925_No19.json",
         "data/processed/GSE285925_log2cpm.csv.gz", "data/processed/GSE285925_sample_metadata.csv",
         "No19", 0, 1, "boundary"),
        ("GSE285925_No20", "results/signature_validation_GSE285925_No20.json",
         "data/processed/GSE285925_log2cpm.csv.gz", "data/processed/GSE285925_sample_metadata.csv",
         "No20", 0, 1, "boundary"),
    ]

    root = args.output.resolve().parents[1]
    rng = np.random.default_rng(args.seed)
    report: dict[str, object] = {
        "scope": (
            "Effect-matched control gene sets for frozen-signature transfer. The signature was "
            "selected on the absolute GSE42783 red-ripe-minus-mature-green effect, so transfer "
            "statistics are reported against controls sharing that selection variable, source "
            "expression level, and gene length. Signature membership, reference direction and "
            "cohort roles are unchanged."
        ),
        "frozen_signature_genes": len(signature),
        "draws": args.draws,
        "seed": args.seed,
        "covariates": ["absolute GSE42783 source effect", "GSE42783 endpoint mean level",
                       "log10 Ensembl Plants release-63 gene length"],
        "rank_control_definition": "genes ranked 501-1000 by absolute source effect, signature members excluded",
        "rank_control_genes": len(rank_control),
        "primary_control": "rank_501_1000",
        "control_construction_note": (
            "The 464 frozen genes are drawn from the 500 largest absolute source effects, so "
            "the signature occupies the extreme tail of the selection variable and a covariate-"
            "balanced non-signature control set does not exist. The matched draws below are "
            "reported for completeness with their balance diagnostics, but they carry "
            "systematically smaller source effects and therefore understate the control; the "
            "rank 501-1000 set is the closest constructible comparison and is treated as primary."
        ),
        "cohorts": {},
    }
    draw_rows: list[dict[str, object]] = []

    for name, frozen_path, matrix_path, metadata_path, genotype, early, late, role in cohorts:
        if frozen_path is None:
            benchmark = pd.read_csv(
                root / "results/external_validation_GSE78733/GSE78733_frozen_external_benchmark.csv"
            )
            frozen_row = benchmark.loc[benchmark["method"].eq("frozen_signature")].iloc[0]
            frozen_report = {
                "direction_concordance": float(frozen_row["direction_concordance"]),
                "spearman_effect_correlation": float(frozen_row["spearman_effect_rank"]),
            }
        else:
            frozen_report = json.loads((root / frozen_path).read_text(encoding="utf-8"))
        matrix = load_matrix(root / matrix_path)
        metadata = pd.read_csv(root / metadata_path)
        imported = endpoint_effect(matrix, metadata, genotype, early, late)

        observed = transfer_statistics(source_effect, imported, signature)
        if observed is None:
            raise ValueError(f"{name}: too few measured signature genes")

        # Reproduce the frozen statistic before adding anything to it.
        for key, frozen_key in (("direction", "direction_concordance"),
                                ("rho", "spearman_effect_correlation")):
            delta = abs(observed[key] - float(frozen_report[frozen_key]))
            if delta > 5e-3:
                raise ValueError(f"{name}: recomputed {key} differs from frozen by {delta:.4f}")

        measured = [g for g in signature if g in imported.index and g in source_effect.index]
        eligible = [g for g in measured if g in z.index]
        candidates = np.array(sorted(
            (set(z.index) & set(imported.index)) - set(signature)
        ))
        if len(candidates) < len(eligible):
            raise ValueError(f"{name}: not enough covariate-complete background genes")
        pools = matched_pools(eligible, candidates, z)

        directions, rhos, smds = [], [], []
        for draw in range(1, args.draws + 1):
            available = np.ones(len(candidates), dtype=bool)
            picked: list[str] = []
            for pool in pools:
                feasible = pool[available[pool]][:10]
                if not len(feasible):
                    feasible = np.flatnonzero(available)
                choice = int(rng.choice(feasible))
                available[choice] = False
                picked.append(str(candidates[choice]))
            stats = transfer_statistics(source_effect, imported, picked)
            directions.append(stats["direction"])
            rhos.append(stats["rho"])
            smds.append(float((z.loc[eligible].mean() - z.loc[picked].mean()).abs().max()))
            draw_rows.append({"cohort": name, "draw": draw,
                              "direction": stats["direction"], "rho": stats["rho"]})

        directions = np.asarray(directions)
        rhos = np.asarray(rhos)
        rank_stats = transfer_statistics(source_effect, imported, rank_control)

        report["cohorts"][name] = {
            "role": role,
            "signature": {"genes": observed["genes"],
                          "direction_concordance": observed["direction"],
                          "spearman_effect_correlation": observed["rho"]},
            "matched_control": {
                "eligible_signature_genes": len(eligible),
                "covariate_complete_background_genes": int(len(candidates)),
                "direction_median": float(np.median(directions)),
                "direction_upper_95": float(np.quantile(directions, 0.95)),
                "direction_empirical_p": float((1 + np.sum(directions >= observed["direction"])) / (1 + args.draws)),
                "rho_median": float(np.median(rhos)),
                "rho_upper_95": float(np.quantile(rhos, 0.95)),
                "rho_empirical_p": float((1 + np.sum(rhos >= observed["rho"])) / (1 + args.draws)),
                "max_absolute_standardized_mean_difference": {
                    "median": float(np.median(smds)),
                    "upper_95_percentile": float(np.quantile(smds, 0.95)),
                    "maximum": float(np.max(smds)),
                },
                "balance_achieved": bool(np.median(smds) <= 0.25),
                "interpretation": (
                    "Balanced: controls are comparable on the selection variable."
                    if np.median(smds) <= 0.25 else
                    "NOT balanced. The signature is the extreme tail of the absolute source "
                    "effect distribution, so no non-signature gene has a comparable effect "
                    "magnitude and the draws carry systematically smaller source effects. "
                    "These draws therefore do not isolate the selection variable and are "
                    "reported only as a lower bound; the rank 501-1000 control is the "
                    "closest constructible comparison."
                ),
            },
            "rank_501_1000_control": (
                {"genes": rank_stats["genes"],
                 "direction_concordance": rank_stats["direction"],
                 "spearman_effect_correlation": rank_stats["rho"]}
                if rank_stats else None
            ),
        }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.null_output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame(draw_rows).to_csv(args.null_output, index=False, compression="infer")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
