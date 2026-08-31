"""Post-review, selection-aware sensitivity analyses over frozen processed inputs.

See docs/hr_revision_20260831/analysis_protocol.md. No source/result overwrite,
network publication, signature reselection, or raw-experiment generation.
"""
from __future__ import annotations
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from matched_transfer_controls import load_matrix, endpoint_effect
from matched_crosslayer_enrichment import truth, ids_from_first_column, core

OUT = ROOT / "results/hr_revision_20260831"
COHORTS = [
    ("GSE267238", "log2cpm", "genotype", "WT", 0, 1, "external_used_in_frozen_grade"),
    ("GSE128739", "log2cpm", "genotype", "Moneymaker", 0, 3, "external_used_in_frozen_grade"),
    ("GSE235023", "log2cpm", "genotype", "AC_exocarp", 0, 1, "external_not_used_in_grade"),
    ("GSE78733", "rank01", "genotype", "WT", 1, 3, "external_not_used_in_grade"),
    ("GSE285925", "log2cpm", "genotype", "No19", 0, 1, "boundary_not_used_in_grade"),
    ("GSE285925", "log2cpm", "genotype", "No20", 0, 1, "boundary_not_used_in_grade"),
    ("GSE183836", "log2cpm", "treatment", "healthy", 0, 3, "development_reuse_sensitivity"),
]

def freeze_baseline():
    target = ROOT / "output/hr_revision_20260831_baseline"
    if target.exists():
        return
    files = [*ROOT.glob("docs/*v1.md"), ROOT / "docs/hr_cover_letter.md",
             ROOT / "docs/hr_review_20260831/synthesis.md",
             ROOT / "config/hr_submission_declarations.json",
             ROOT / "results/hr_submission/rin_candidate_evidence_table.csv",
             ROOT / "results/ripening_atlas_resource_v3.csv"]
    for source in files:
        dest = target / source.relative_to(ROOT)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
    manifest = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(target)}"
                for p in sorted(target.rglob("*")) if p.is_file()]
    (target / "SHA256SUMS.txt").write_text("\n".join(manifest) + "\n")

def bootstrap_effects(early: np.ndarray, late: np.ndarray) -> np.ndarray:
    """All ordered endpoint column resamples, each equally probable.

    Whole biological columns are sampled together. Rows are genes. The input
    transformations are held fixed, and source selection is not re-estimated.
    """
    early_means = np.stack([early[:, idx].mean(axis=1) for idx in
                           itertools.product(range(early.shape[1]), repeat=early.shape[1])], axis=1)
    late_means = np.stack([late[:, idx].mean(axis=1) for idx in
                          itertools.product(range(late.shape[1]), repeat=late.shape[1])], axis=1)
    return (late_means[:, None, :] - early_means[:, :, None]).reshape(early.shape[0], -1)

def loo_effects(early: np.ndarray, late: np.ndarray) -> np.ndarray:
    return np.stack([np.delete(late, j, axis=1).mean(axis=1) -
                     np.delete(early, i, axis=1).mean(axis=1)
                     for i in range(early.shape[1]) for j in range(late.shape[1])], axis=1)

def rank_corr(reference: np.ndarray, effects: np.ndarray) -> np.ndarray:
    a = pd.Series(reference).rank().to_numpy()
    b = pd.DataFrame(effects).rank(axis=0).to_numpy()
    a = a - a.mean()
    b = b - b.mean(axis=0)
    return (a[:, None] * b).sum(axis=0) / np.sqrt((a*a).sum() * (b*b).sum(axis=0))

def bh(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    adjusted = np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values)+1))[::-1])[::-1]
    result = np.empty(len(values))
    result[order] = np.minimum(adjusted, 1)
    return result

def source_inputs():
    effects = []
    for accession in ("GSE42783", "GSE108415"):
        matrix = load_matrix(ROOT / f"data/processed/{accession}_rank01.csv.gz")
        metadata = pd.read_csv(ROOT / f"data/processed/{accession}_sample_metadata.csv")
        effects.append(endpoint_effect(matrix, metadata, "WT", 0, 3))
    sources = pd.concat(effects, axis=1, keys=["source_a", "source_b"]).dropna()
    sources["source_b_confirms"] = sources.source_a * sources.source_b > 0
    return sources

def selection_support(sources, signature):
    protein = pd.read_csv(ROOT / "results/PXD051570_gene_evidence.csv")
    compact_universe = ROOT / "data/processed/GPL15968_gene_universe.csv"
    occupancy = (set(pd.read_csv(compact_universe).gene_id) if compact_universe.exists()
                 else ids_from_first_column(ROOT / "data/raw/GPL15968_110707_RDKK310_Slyc_ChIP.ndf.gz"))
    binding = pd.read_csv(ROOT / "results/GSE40257_rin_targets.csv")
    admitted = set(binding.loc[truth(binding.admitted_binding), "gene_id"]) & occupancy
    de = pd.read_csv(ROOT / "results/GSE210589/deseq2/GSE210589_rin1_vs_WT_DESeq2.csv")
    universes = {
        "protein_differential": set(protein.loc[truth(protein.protein_measured), "gene_id"]),
        "RIN_occupancy": occupancy,
        "rin1_differential_among_occupied_testable_genes": admitted & set(de.gene_id),
    }
    rows = []
    for name, universe in universes.items():
        observed = sorted(signature & universe & set(sources.index))
        background_all = sorted((universe & set(sources.index)) - signature)
        background = [g for g in background_all if bool(sources.loc[g, "source_b_confirms"])]
        a = sources.loc[observed, "source_a"].abs()
        b = sources.loc[background, "source_a"].abs()
        scale = pd.concat([a, b]).std(ddof=1)
        best_mean = b.nlargest(len(a)).mean() if len(b) >= len(a) else np.nan
        rows.append({"layer": name, "assay_universe_n": len(universe),
                     "signature_measured_n": len(signature & universe),
                     "signature_with_both_sources_n": len(a),
                     "background_with_both_sources_n": len(background_all),
                     "background_source_b_confirmed_n": len(b),
                     "signature_absolute_source_a_min": a.min(),
                     "signature_absolute_source_a_mean": a.mean(),
                     "background_absolute_source_a_max": b.max(),
                     "best_equal_size_background_absolute_source_a_mean": best_mean,
                     "absolute_source_a_standardized_mean_gap_lower_bound": (a.mean()-best_mean)/scale,
                     "signature_in_background_empirical_range_n": int(a.between(b.min(), b.max()).sum()),
                     "range_overlap": bool(a.min() <= b.max() and b.min() <= a.max()),
                     "selection_adjusted_enrichment_estimable": False,
                     "interpretation": "No common empirical source-effect support for the full frozen signature; no selection-adjusted P value is estimated."})
        if rows[-1]["range_overlap"]:
            raise RuntimeError("Unexpected empirical overlap; inspect selection support before inference")
    pd.DataFrame(rows).to_csv(OUT / "selection_support.csv", index=False)
    old = json.loads((ROOT / "results/matched_crosslayer_enrichment.json").read_text())
    raw = [v["matched_null"]["matched_null_empirical_p"] for v in old["layers"].values()]
    table = [{"layer": name, "nominal_p": p, "exploratory_bh_three_layers": q,
              "selection_adjusted": False,
              "matching_variables": "GSE267238 detection, absolute endpoint effect, gene length"}
             for (name, _), p, q in zip(old["layers"].items(), raw, bh(raw))]
    pd.DataFrame(table).to_csv(OUT / "crosslayer_multiplicity.csv", index=False)

def main():
    freeze_baseline()
    OUT.mkdir(parents=True, exist_ok=True)
    sources = source_inputs()
    signature = set(pd.read_csv(ROOT / "results/concordant_mature_green_to_ripe_signature.csv").gene_id)
    ordered = sources.source_a.abs().sort_values(ascending=False, kind="stable").index.tolist()
    derived = {g for g in ordered[:500] if sources.loc[g, "source_b_confirms"]}
    assert derived == signature and len(signature) == 464
    methods = {"frozen_signature": sorted(signature), "source_top500": ordered[:500],
               "source_top464": ordered[:464], "rank501_1000": ordered[500:1000],
               "rank501_1000_confirmed": [g for g in ordered[500:1000] if sources.loc[g, "source_b_confirms"]]}
    pd.DataFrame([{"method": method, "gene_id": gene} for method, genes in methods.items()
                  for gene in genes]).to_csv(OUT / "fixed_comparison_sets.csv", index=False)
    candidates = pd.read_csv(ROOT / "results/hr_submission/rin_candidate_evidence_table.csv")
    transfer_rows, candidate_rows, sample_rows, bootstrap_rows = [], [], [], []
    for accession, scale, field, value, early_stage, late_stage, role in COHORTS:
        context = accession + ("_"+value if accession in {"GSE285925", "GSE183836"} else "")
        matrix = load_matrix(ROOT / f"data/processed/{accession}_{scale}.csv.gz")
        metadata = pd.read_csv(ROOT / f"data/processed/{accession}_sample_metadata.csv")
        selected = metadata[metadata[field].eq(value)]
        early_ids = selected.loc[selected.stage_ordinal.eq(early_stage), "sample_id"].tolist()
        late_ids = selected.loc[selected.stage_ordinal.eq(late_stage), "sample_id"].tolist()
        assert len(early_ids) >= 2 and len(late_ids) >= 2
        for stage, ids in (("early", early_ids), ("late", late_ids)):
            for sample in ids:
                sample_rows.append({"context": context, "role": role, "endpoint": stage,
                                    "sample_id": sample, "scale": scale,
                                    "unit": "biological sample; technical lanes already collapsed where applicable"})
        genes = sorted(set(sources.index) & set(matrix.index))
        matrix = matrix.loc[genes]
        source = sources.loc[genes, "source_a"].to_numpy()
        early = matrix[early_ids].to_numpy()
        late = matrix[late_ids].to_numpy()
        point = late.mean(axis=1) - early.mean(axis=1)
        boot = bootstrap_effects(early, late)
        loo = loo_effects(early, late)
        idx = {g:i for i,g in enumerate(genes)}
        method_stats = {}
        for method, members in methods.items():
            take = np.array([idx[g] for g in members if g in idx], dtype=int)
            sign = np.sign(source[take])
            directions = (sign[:,None] == np.sign(boot[take])).mean(axis=0)
            correlations = rank_corr(source[take], boot[take])
            loo_directions = (sign[:,None] == np.sign(loo[take])).mean(axis=0)
            method_stats[method] = (directions, correlations)
            record = {"context": context, "role": role, "method": method,
                      "n_early": len(early_ids), "n_late": len(late_ids), "measured_genes": len(take),
                      "direction": float((sign == np.sign(point[take])).mean()),
                      "rho": float(rank_corr(source[take], point[take,None])[0]),
                      "bootstrap_ordered_resamples": boot.shape[1],
                      "direction_bootstrap_p025": np.quantile(directions,.025),
                      "direction_bootstrap_p975": np.quantile(directions,.975),
                      "rho_bootstrap_p025": np.quantile(correlations,.025),
                      "rho_bootstrap_p975": np.quantile(correlations,.975),
                      "leave_one_each_endpoint_min": loo_directions.min(),
                      "leave_one_each_endpoint_max": loo_directions.max()}
            transfer_rows.append(record)
        for method, (direction, rho) in method_stats.items():
            delta = method_stats["frozen_signature"][0] - direction
            record = next(r for r in transfer_rows if r["context"]==context and r["method"]==method)
            record.update(signature_minus_method_direction_bootstrap_p025=np.quantile(delta,.025),
                          signature_minus_method_direction_bootstrap_p975=np.quantile(delta,.975))
            for k, (d,r,dd) in enumerate(zip(direction,rho,delta)):
                bootstrap_rows.append({"context":context,"method":method,"resample":k+1,
                                       "direction":d,"rho":r,"signature_minus_method_direction":dd})
        for _, row in candidates.iterrows():
            gene = row.gene_id
            record = {"gene_id":gene,"display_name":row.display_name,"frozen_grade":row.candidate_grade,
                      "context":context,"role":role,"n_early":len(early_ids),"n_late":len(late_ids),
                      "effect_scale":scale,"measured":gene in idx}
            if gene in idx:
                j = idx[gene]
                record.update(effect=point[j], reference_effect=source[j],
                              same_direction=bool(np.sign(point[j])==np.sign(source[j])),
                              bootstrap_same_direction_frequency=float((np.sign(boot[j])==np.sign(source[j])).mean()),
                              effect_bootstrap_p025=np.quantile(boot[j],.025),
                              effect_bootstrap_p975=np.quantile(boot[j],.975))
            candidate_rows.append(record)
        print(context, "complete", flush=True)
    pd.DataFrame(transfer_rows).to_csv(OUT / "transfer_comparisons.csv",index=False)
    pd.DataFrame(candidate_rows).to_csv(OUT / "candidate_external_contexts.csv",index=False)
    pd.DataFrame(sample_rows).to_csv(OUT / "resampling_units.csv",index=False)
    pd.DataFrame(bootstrap_rows).to_csv(OUT / "bootstrap_statistics.csv.gz",index=False,compression="gzip")
    selection_support(sources, signature)
    summary = {"frozen_signature":464,"comparison_set_sizes":{k:len(v) for k,v in methods.items()},
               "candidate_rows":len(candidate_rows),"external_studies":4,"boundary_studies":1,
               "development_reuse_sensitivity_studies":1,
               "resampling":"exact empirical endpoint-column bootstrap conditional on processed inputs and frozen reference; not a raw-pipeline or across-study confidence interval",
               "source_selection_rerun":False,"raw_data_generated":False,
               "protocol_sha256":hashlib.sha256((ROOT / "docs/hr_revision_20260831/analysis_protocol.md").read_bytes()).hexdigest()}
    (OUT / "analysis_summary.json").write_text(json.dumps(summary,indent=2)+"\n")

if __name__ == "__main__":
    main()
