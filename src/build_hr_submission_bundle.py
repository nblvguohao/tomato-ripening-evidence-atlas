#!/usr/bin/env python3
"""Build a resource-led Horticulture Research submission bundle.

This script does not refit models or alter frozen atlas v3 results. It derives a
submission-facing candidate table, a compact summary, and six SVG figures from
already frozen outputs. StudyShield is intentionally excluded from candidate
grading; its benchmark remains a supplementary negative result.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COLORS = {
    # Restrained, colour-blind-aware palette with neutral typography. The
    # palette is intentionally shared by all six figures to approximate the
    # compact, low-decoration visual system used in Nature-style line art.
    "navy": "#202124",
    "blue": "#3C78A8",
    "teal": "#008F7A",
    "green": "#4C956C",
    "orange": "#D56A1A",
    "red": "#BD4B43",
    "purple": "#7356A5",
    "gray": "#6E7378",
    "grid": "#D9DDE1",
    "light": "#F1F3F5",
    "pale": "#FFFFFF",
}

PENDING_LITERATURE_STATUSES = {
    "not_screened",
    "pending_full_alias_and_full_text_screen",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def as_bool(value: object) -> bool:
    return str(value).strip().lower() == "true"


def as_float(value: object, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def display_name(row: dict[str, str]) -> str:
    return row.get("ensembl_gene_name", "").strip() or row["gene_id"]


def annotation_theme(row: dict[str, str]) -> tuple[str, str]:
    """Assign a descriptive, annotation-derived theme without enrichment claims."""
    gene = row["gene_id"]
    name = row.get("ensembl_gene_name", "").strip()
    description = row.get("gene_description", "").lower()
    text = f"{name.lower()} {description}"
    specific = {
        "Solyc01g095080": ("Ethylene and ripening transcription", "ACS2 annotation"),
        "Solyc10g006880": ("Ethylene and ripening transcription", "NOR annotation"),
        "Solyc08g080630": ("Ethylene and ripening transcription", "ethylene-responsive annotation"),
        "Solyc08g005610": ("Hormone and growth signalling", "CYP707A2 annotation"),
        "Solyc08g005770": ("Flavour and specialised metabolism", "AAT annotation"),
        "Solyc01g006540": ("Flavour and specialised metabolism", "loxC annotation"),
        "Solyc08g014130": ("Flavour and specialised metabolism", "IPMS2 annotation"),
        "Solyc10g084600": ("Plastid light adaptation", "published SlBGH2B functional study; see grade-A literature screen"),
    }
    if gene in specific:
        return specific[gene]
    rules = [
        ("Cell wall and proteolysis", ("expansin", "pectin methylesterase", "carboxypeptidase")),
        ("Hormone and growth signalling", ("auxin", "argos", "flagellin sensing", "receptor-like", "bhlh transcription")),
        ("Flavour and specialised metabolism", ("cytochrome p450", "oxidoreductase", "aldo-keto")),
        ("Defence and stress response", ("pathogen-related", "wound-responsive", "selenium-binding", "lipid-transfer", "proteinase inhibitor")),
        ("Photosynthetic and organellar metabolism", ("chlorophyll a-b", "nudix", "iron-sulfur", "nfu1")),
        ("Genome and protein homeostasis", ("deoxyuridine", "poly [adp-ribose]", "dna polymerase", "wd40", "ubiquitin", "big1")),
        ("Limited functional annotation", ("uncharacterized", "hypothetical", "unknown protein", "duf")),
    ]
    for theme, keywords in rules:
        for keyword in keywords:
            if keyword in text:
                return theme, f"matched annotation keyword: {keyword}"
    if not name and not description.strip():
        return "Limited functional annotation", "no functional description in frozen annotation"
    return "Other annotated functions", "annotation does not match a prespecified descriptive theme"


def candidate_grade(row: dict[str, str]) -> tuple[str, str]:
    consensus = as_bool(row.get("all_four_same_direction"))
    protein = row.get("protein_support_status", "")
    protein_concordant = protein in {
        "differential_direction_concordant",
        "measured_direction_concordant",
    }
    if not consensus:
        return (
            "C",
            "RIN binding plus rin-1 response, but incomplete four-cohort expression replication",
        )
    if protein_concordant:
        return (
            "A",
            "four-cohort expression consensus plus RIN binding, rin-1 response, and concordant protein direction",
        )
    return (
        "B",
        "four-cohort expression consensus plus RIN binding and rin-1 response; no concordant protein support",
    )


def literature_screen_complete(status: object) -> bool:
    """Return whether a literature-screen row contains a completed judgement.

    A registry row created for auditability is not itself a completed literature
    review. Keeping this distinction explicit prevents pending rows from being
    reported as evidence in the submission summary.
    """
    return str(status).strip() not in PENDING_LITERATURE_STATUSES


def build_candidates(
    atlas_rows: list[dict[str, str]],
    perturbation_rows: list[dict[str, str]],
    literature_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    perturbation = {row["gene_id"]: row for row in perturbation_rows}
    literature = {row["gene_id"]: row for row in literature_rows}
    candidates: list[dict[str, object]] = []
    for atlas in atlas_rows:
        record = perturbation.get(atlas["gene_id"])
        if not record or record.get("evidence_grade") != "binding_plus_independent_perturbation_expression":
            continue
        grade, definition = candidate_grade(atlas)
        theme, theme_basis = annotation_theme(atlas)
        protein_status = atlas.get("protein_support_status", "not_measured")
        protein_concordant = protein_status in {
            "differential_direction_concordant",
            "measured_direction_concordant",
        }
        description = atlas.get("gene_description", "")
        limited_functional_annotation = theme == "Limited functional annotation"
        anchor = atlas["gene_id"] in {
            "Solyc01g095080",  # ACS2
            "Solyc10g006880",  # NOR
            "Solyc06g051800",  # EXP1
            "Solyc08g005770",  # AAT
            "Solyc01g006540",  # loxC
            "Solyc08g005610",  # CYP707A2
        }
        tissue_count = int(as_float(atlas.get("tissue_context_count"), 0))
        tissue_consistency = as_float(atlas.get("tissue_context_direction_consistency"), 0)
        literature_record = literature.get(atlas["gene_id"], {})
        row: dict[str, object] = {
            "gene_id": atlas["gene_id"],
            "display_name": display_name(atlas),
            "candidate_grade": grade,
            "grade_definition": definition,
            "annotation_theme": theme,
            "theme_basis": theme_basis,
            "biological_anchor": anchor,
            "limited_functional_annotation": limited_functional_annotation,
            "four_cohort_direction_consensus": as_bool(atlas.get("all_four_same_direction")),
            "ripening_direction": "higher" if as_float(atlas.get("effect_source_a"), 0) > 0 else "lower",
            "rin1_response": record.get("perturbation_direction", ""),
            "rin1_vs_wt_log2_fold_change": record.get("log2FoldChange", ""),
            "rin1_vs_wt_bh_fdr": record.get("padj", ""),
            "rin_binding_replicates": record.get("binding_replicates", ""),
            "protein_support_status": protein_status,
            "protein_direction_concordant": protein_concordant,
            "protein_differential": as_bool(atlas.get("protein_differential_BR8_vs_MG")),
            "phosphosite_differential": as_bool(atlas.get("phosphosite_differential_BR8_vs_MG")),
            "tissue_context_count": tissue_count,
            "tissue_direction_consistency": tissue_consistency,
            "gene_description": description,
            "claim_boundary": "candidate priority only; no direct causality or regulatory sign inferred",
            "literature_screen_status": literature_record.get("literature_screen_status", "not_screened"),
            "literature_screen_scope": literature_record.get("literature_screen_scope", ""),
            "literature_screen_source": literature_record.get("literature_screen_source", ""),
            "literature_screen_url": literature_record.get("literature_screen_url", ""),
            "literature_screen_conclusion": literature_record.get("literature_screen_conclusion", ""),
            "literature_screen_date": literature_record.get("screen_date", ""),
        }
        candidates.append(row)
    order = {"A": 0, "B": 1, "C": 2}
    candidates.sort(
        key=lambda row: (
            order[str(row["candidate_grade"])],
            not bool(row["biological_anchor"]),
            not bool(row["limited_functional_annotation"]),
            not bool(row["protein_differential"]),
            not bool(row["phosphosite_differential"]),
            -abs(as_float(row["rin1_vs_wt_log2_fold_change"], 0)),
            str(row["gene_id"]),
        )
    )
    for rank, row in enumerate(candidates, start=1):
        row["submission_priority_rank"] = rank
    return candidates


def esc(value: object) -> str:
    return html.escape(str(value))


def svg_text(x: float, y: float, value: object, size: int = 14, anchor: str = "start", color: str = "#202124", weight: str = "normal") -> str:
    return f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" font-weight="{weight}" fill="{color}">{esc(value)}</text>'


def svg_begin(title: str, subtitle: str, width: int = 1100, height: int = 720) -> list[str]:
    output_height = height - 75
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{output_height}" viewBox="0 0 {width} {output_height}">',
        f'<title>{html.escape(title)}</title>',
        f'<desc>{html.escape(subtitle)}</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<g transform="translate(0,-75)">',
    ]


def svg_write(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines + ["</g>", "</svg>", ""]), encoding="utf-8")


def figure1_design(summary: dict[str, object], path: Path) -> None:
    lines = svg_begin(
        "Figure 1. Frozen discovery-to-evidence design",
        "Study-level independence is preserved; cross-modality evidence annotates rather than reselects the 464-gene signature.",
    )
    boxes = [
        (70, 175, 190, 105, "2 source studies", "frozen 464-gene\nMG-to-ripe signature", COLORS["blue"]),
        (325, 175, 190, 105, "independent cohorts", "cross-platform\nexpression transfer", COLORS["teal"]),
        (580, 175, 190, 105, "evidence layers", "protein, RIN binding,\nrin-1 perturbation", COLORS["orange"]),
        (835, 175, 190, 105, "released resource", "gene-level evidence\nand claim boundaries", COLORS["purple"]),
    ]
    lines.append(svg_text(45, 120, "a", 18, weight="bold"))
    for index, (x, y, width, height, heading, body, color) in enumerate(boxes):
        lines.append(f'<rect x="{x}" y="{y}" width="{width}" height="{height}" fill="white" stroke="{COLORS["grid"]}" stroke-width="1.5"/>')
        lines.append(f'<rect x="{x}" y="{y}" width="{width}" height="5" fill="{color}"/>')
        lines.append(svg_text(x + width / 2, y + 32, heading, 15, anchor="middle", weight="bold"))
        for line_index, body_line in enumerate(body.split("\n")):
            lines.append(svg_text(x + width / 2, y + 64 + line_index * 21, body_line, 13, anchor="middle"))
        if index < len(boxes) - 1:
            lines.append(f'<line x1="{x+width}" y1="{y+height/2}" x2="{x+width+58}" y2="{y+height/2}" stroke="{COLORS["gray"]}" stroke-width="1.5"/>')
            lines.append(f'<polygon points="{x+width+58},{y+height/2} {x+width+48},{y+height/2-5} {x+width+48},{y+height/2+5}" fill="{COLORS["gray"]}"/>')
    metrics = [
        ("464", "frozen genes"),
        ("33", "four-cohort consensus + RIN response"),
        (str(summary["grade_counts"]["A"]), "grade A candidates"),
        ("6", "main figures"),
    ]
    lines.append(f'<line x1="70" y1="405" x2="1030" y2="405" stroke="{COLORS["grid"]}" stroke-width="1"/>')
    for index, (value, label) in enumerate(metrics):
        x = 140 + index * 255
        lines.append(svg_text(x, 470, value, 31, anchor="middle", weight="bold"))
        lines.append(svg_text(x, 502, label, 12, anchor="middle", color=COLORS["gray"]))
    lines.append(f'<line x1="45" y1="575" x2="1055" y2="575" stroke="{COLORS["grid"]}" stroke-width="1"/>')
    lines.append(svg_text(45, 615, "Every added layer annotates the frozen signature; none changes membership or reference direction.", 14, weight="bold"))
    svg_write(path, lines)


def figure2_transfer(root: Path, external_rows: list[dict[str, str]], path: Path) -> None:
    validations = [
        ("GSE267238 pericarp", "signature_validation_GSE267238.json"),
        ("GSE128739 time course", "signature_validation_GSE128739.json"),
        ("GSE235023 exocarp", "signature_validation_GSE235023.json"),
        ("GSE285925 No19", "signature_validation_GSE285925_No19.json"),
        ("GSE285925 No20", "signature_validation_GSE285925_No20.json"),
    ]
    data: list[tuple[str, float, float]] = []
    for label, filename in validations:
        record = json.loads((root / "results" / filename).read_text())
        data.append((label, float(record["spearman_effect_correlation"]), float(record["direction_concordance"])))
    external = next(row for row in external_rows if row["method"] == "frozen_signature")
    data.insert(3, ("GSE78733 untouched", as_float(external["spearman_effect_rank"]), as_float(external["direction_concordance"])))
    lines = svg_begin(
        "Figure 2. Frozen expression signature transfers across independent cohorts",
        "Effect-rank correlation and direction agreement are shown for a signature imported without gene reselection.",
        height=790,
    )
    lines.append(svg_text(45, 112, "a", 18, weight="bold"))
    lines.append(svg_text(310, 112, "Effect-rank Spearman", 13, weight="bold"))
    lines.append(svg_text(790, 112, "Direction agreement", 13, weight="bold"))
    for index, (label, rho, direction) in enumerate(data):
        y = 150 + index * 92
        boundary = "285925" in label
        rho_color = COLORS["orange"] if boundary else COLORS["blue"]
        direction_color = COLORS["orange"] if boundary else COLORS["teal"]
        lines.append(svg_text(55, y + 17, label, 14, weight="bold" if boundary else "normal"))
        lines.append(f'<rect x="310" y="{y}" width="430" height="22" fill="{COLORS["light"]}"/>')
        lines.append(f'<rect x="310" y="{y}" width="{max(rho, 0)*430:.1f}" height="22" fill="{rho_color}"/>')
        lines.append(svg_text(755, y + 17, f"{rho:.3f}", 13, anchor="end", weight="bold"))
        lines.append(f'<rect x="790" y="{y}" width="240" height="22" fill="{COLORS["light"]}"/>')
        lines.append(f'<rect x="790" y="{y}" width="{direction*240:.1f}" height="22" fill="{direction_color}"/>')
        lines.append(svg_text(1050, y + 17, f"{direction:.1%}", 13, anchor="end", weight="bold"))
    lines.append(f'<line x1="45" y1="690" x2="1055" y2="690" stroke="{COLORS["grid"]}" stroke-width="1"/>')
    lines.append(svg_text(45, 728, "The prospective breaker-to-mature GSE285925 contrasts define the boundary of transfer and were retained without retuning.", 14, weight="bold"))
    svg_write(path, lines)


def figure3_orthogonal(protein: dict[str, object], binding: dict[str, object], perturbation: dict[str, object], path: Path) -> None:
    lines = svg_begin(
        "Figure 3. Cross-modality protein and RIN evidence supports the frozen signature",
        "Each layer is evaluated against its own measured background; correspondence is not labelled causal.",
    )
    panels = [
        (65, 135, "Protein abundance", COLORS["green"], [
            (f"{protein['signature_protein_genes_measured']}/464", "signature proteins measured"),
            (f"{float(protein['signature_protein_direction_concordance']):.1%}", "transcript-protein direction agreement"),
            (f"OR {float(protein['protein_differential_enrichment']['odds_ratio']):.2f}", "differential-protein enrichment"),
        ]),
        (395, 135, "RIN occupancy", COLORS["orange"], [
            (f"{binding['signature_enrichment']['signature_bound_genes']}/{binding['signature_enrichment']['signature_genes_in_array_universe']}", "signature genes bound"),
            (f"OR {float(binding['signature_enrichment']['odds_ratio']):.2f}", "binding enrichment"),
            (f"p={float(binding['signature_enrichment']['fisher_p_value']):.1e}", "measured-universe Fisher test"),
        ]),
        (725, 135, "Independent rin-1 response", COLORS["purple"], [
            (f"{perturbation['signature_RIN_bound_differential']}/{perturbation['signature_RIN_bound_differential'] + perturbation['signature_RIN_bound_not_differential']}", "testable bound-signature responses"),
            (f"OR {float(perturbation['signature_vs_other_RIN_target_differential_odds_ratio']):.2f}", "versus other RIN-bound genes"),
            (f"{perturbation['high_grade_lower_in_rin1']} lower / {perturbation['high_grade_higher_in_rin1']} higher", "direction in rin-1"),
        ]),
    ]
    for panel_index, (x, y, heading, color, metrics) in enumerate(panels):
        lines.append(f'<rect x="{x}" y="{y}" width="285" height="390" fill="white" stroke="{COLORS["grid"]}" stroke-width="1.5"/>')
        lines.append(f'<rect x="{x}" y="{y}" width="285" height="5" fill="{color}"/>')
        lines.append(svg_text(x + 18, y + 35, chr(ord("a") + panel_index), 17, weight="bold"))
        lines.append(svg_text(x + 48, y + 35, heading, 16, weight="bold"))
        for index, (value, label) in enumerate(metrics):
            yy = y + 112 + index * 105
            lines.append(svg_text(x + 142, yy, value, 25, anchor="middle", weight="bold"))
            lines.append(svg_text(x + 142, yy + 28, label, 12, anchor="middle", color=COLORS["gray"]))
    phospho_p = float(protein["phosphosite_differential_enrichment"]["matched_null_empirical_p"])
    lines.append(f'<line x1="65" y1="575" x2="1030" y2="575" stroke="{COLORS["grid"]}" stroke-width="1"/>')
    lines.append(svg_text(65, 610, f"Phosphosite changes are retained gene by gene but are not used as an enrichment claim (matched-null p={phospho_p:.3f}).", 14, color=COLORS["gray"]))
    lines.append(svg_text(65, 645, "Binding plus perturbation-expression change prioritises candidates; it does not establish direct regulation, sign, or mechanism.", 14, weight="bold"))
    svg_write(path, lines)


def figure4_applicability(tissue: dict[str, object], root: Path, path: Path) -> None:
    contexts = sorted(tissue["contexts"], key=lambda row: float(row["direction_concordance"]), reverse=True)
    lines = svg_begin(
        "Figure 4. Tissue applicability is broad but not uniform",
        "SRP109982 supplies replicate-averaged M82 profiles for applicability mapping, not model fitting.",
        height=850,
    )
    lines.append(svg_text(45, 100, "a", 18, weight="bold"))
    lines.append(svg_text(250, 100, "Direction agreement", 13, weight="bold"))
    lines.append(svg_text(845, 100, "Effect-rank correlation", 13, weight="bold"))
    for index, record in enumerate(contexts):
        y = 120 + index * 57
        direction = float(record["direction_concordance"])
        rho = float(record["spearman_effect_correlation"])
        lines.append(svg_text(55, y + 16, record["context"], 13))
        lines.append(f'<line x1="250" y1="{y+11}" x2="750" y2="{y+11}" stroke="{COLORS["grid"]}" stroke-width="1"/>')
        lines.append(f'<rect x="250" y="{y}" width="{direction*500:.1f}" height="22" fill="{COLORS["teal"]}"/>')
        lines.append(svg_text(770, y + 16, f"{direction:.1%}", 12, anchor="end", weight="bold"))
        lines.append(f'<circle cx="{845 + rho*180:.1f}" cy="{y+11}" r="7" fill="{COLORS["orange"]}"/>')
        lines.append(svg_text(1050, y + 16, f"rho={rho:.3f}", 12, anchor="end"))
    lines.append(svg_text(250, 785, "direction agreement", 12, color=COLORS["gray"]))
    lines.append(svg_text(845, 785, "effect-rank correlation", 12, color=COLORS["gray"]))
    lines.append(svg_text(55, 820, "Individual libraries were unavailable in the processed supplement; uncertainty is therefore not estimated at the replicate level.", 13, color=COLORS["gray"]))
    svg_write(path, lines)


def figure5_modules(candidates: list[dict[str, object]], go_rows: list[dict[str, str]], path: Path) -> None:
    themes = Counter(str(row["annotation_theme"]) for row in candidates)
    significant = [row for row in go_rows if as_float(row.get("bh_fdr"), 1) < .05]
    tested = len(go_rows)
    min_fdr = min((as_float(row.get("bh_fdr"), 1) for row in go_rows), default=math.nan)
    matched_count = int(as_float(go_rows[0].get("matched_signature_gene_count"), 0)) if go_rows else 0
    frozen_count = int(as_float(go_rows[0].get("full_frozen_signature_gene_count"), 0)) if go_rows else 0
    lines = svg_begin(
        "Figure 5. Annotation themes and the matched-null Gene Ontology result",
        "Themes are descriptive labels; Gene Ontology testing uses the exact covariate-eligible signature frame and yields no FDR-significant term.",
        width=1100,
        height=720,
    )
    lines.append(svg_text(45, 92, "a", 18, weight="bold"))
    lines.append(svg_text(78, 92, "Descriptive annotation themes", 16, weight="bold"))
    max_count = max(themes.values())
    for index, (theme, count) in enumerate(sorted(themes.items(), key=lambda item: (-item[1], item[0]))):
        y = 135 + index * 53
        short_theme = theme.replace("Photosynthetic and organellar metabolism", "Photosynthesis / organelles").replace("Ethylene and ripening transcription", "Ethylene / ripening transcription").replace("Genome and protein homeostasis", "Genome / protein homeostasis")
        lines.append(svg_text(45, y + 17, short_theme, 14))
        lines.append(f'<rect x="305" y="{y}" width="220" height="20" fill="{COLORS["light"]}"/>')
        lines.append(f'<rect x="305" y="{y}" width="{count/max_count*220:.1f}" height="20" fill="{COLORS["blue"]}"/>')
        lines.append(svg_text(550, y + 18, count, 14, anchor="end", weight="bold"))
    lines.append(f'<line x1="600" y1="105" x2="600" y2="585" stroke="{COLORS["grid"]}" stroke-width="1"/>')
    lines.append(svg_text(640, 92, "b", 18, weight="bold"))
    lines.append(svg_text(673, 92, "Matched-null GO test", 16, weight="bold"))
    lines.append(svg_text(820, 250, f"{len(significant)}", 58, anchor="middle", weight="bold"))
    lines.append(svg_text(820, 285, f"of {tested} tested terms at FDR < 0.05", 16, anchor="middle"))
    lines.append(svg_text(820, 365, f"minimum FDR = {min_fdr:.3f}", 18, anchor="middle", weight="bold"))
    lines.append(svg_text(820, 420, f"matched frame: {matched_count}/{frozen_count} genes", 16, anchor="middle", color=COLORS["gray"]))
    lines.append(svg_text(820, 455, "1,000 covariate-matched draws", 16, anchor="middle", color=COLORS["gray"]))
    lines.append(f'<line x1="45" y1="590" x2="1055" y2="590" stroke="{COLORS["grid"]}" stroke-width="1"/>')
    lines.append(svg_text(45, 635, "Theme counts are descriptive, not enrichment tests. No signature-wide Gene Ontology enrichment is claimed.", 14, weight="bold"))
    svg_write(path, lines)


def figure6_candidates(candidates: list[dict[str, object]], path: Path) -> None:
    grade_counts = Counter(str(row["candidate_grade"]) for row in candidates)
    four_cohort = sum(bool(row["four_cohort_direction_consensus"]) for row in candidates)
    protein = sum(bool(row["protein_direction_concordant"]) for row in candidates)
    broad_tissue = sum(int(row["tissue_context_count"]) >= 9 and float(row["tissue_direction_consistency"]) >= .8 for row in candidates)
    anchors = [row for row in candidates if bool(row["biological_anchor"])]
    lines = svg_begin(
        "Figure 6. Evidence-tiered summary of 35 RIN-associated candidates",
        "Aggregate grades and evidence counts summarise the candidate set; the complete 35-row matrix is provided in Supplementary Table S6.",
        width=1100,
        height=720,
    )
    grade_colors = {"A": COLORS["green"], "B": COLORS["orange"], "C": COLORS["red"]}
    lines.append(svg_text(45, 92, "a", 18, weight="bold"))
    lines.append(svg_text(78, 92, "Evidence tiers", 16, weight="bold"))
    for index, grade in enumerate(("A", "B", "C")):
        x = 55 + index * 340
        lines.append(f'<rect x="{x}" y="135" width="305" height="170" fill="white" stroke="{COLORS["grid"]}" stroke-width="1.5"/>')
        lines.append(f'<rect x="{x}" y="135" width="305" height="5" fill="{grade_colors[grade]}"/>')
        lines.append(svg_text(x + 152, 190, f"Grade {grade}", 19, anchor="middle", weight="bold"))
        lines.append(svg_text(x + 152, 250, grade_counts[grade], 52, anchor="middle", weight="bold"))
        definition = {"A": "4-cohort + protein", "B": "4-cohort; no protein support", "C": "partial expression replication"}[grade]
        lines.append(svg_text(x + 152, 282, definition, 14, anchor="middle", color=COLORS["gray"]))
    metrics = [
        ("35", "RIN occupancy + rin-1 response"),
        (str(four_cohort), "four-cohort direction consensus"),
        (str(protein), "concordant protein direction"),
        (str(broad_tissue), "broad tissue applicability"),
    ]
    for index, (value, label) in enumerate(metrics):
        x = 55 + index * 260
        lines.append(svg_text(x + 105, 385, value, 34, anchor="middle", weight="bold"))
        lines.append(svg_text(x + 105, 417, label, 13, anchor="middle", color=COLORS["gray"]))
    anchor_labels = ", ".join(str(row["display_name"]) for row in anchors)
    lines.append(f'<rect x="55" y="470" width="985" height="105" fill="{COLORS["light"]}"/>')
    lines.append(svg_text(80, 512, "Biological calibration anchors", 17, color=COLORS["navy"], weight="bold"))
    lines.append(svg_text(80, 548, anchor_labels, 16, color=COLORS["navy"]))
    lines.append(svg_text(55, 625, "Grades prioritise follow-up; they do not establish direct regulation, regulatory sign, or mechanism. Full evidence: Table S6.", 14, weight="bold"))
    svg_write(path, lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas", type=Path, default=ROOT / "results/ripening_atlas_resource_v3.csv")
    parser.add_argument("--perturbation", type=Path, default=ROOT / "results/GSE210589/deseq2/GSE210589_RIN_binding_perturbation_evidence.csv")
    parser.add_argument("--literature-screen", type=Path, default=ROOT / "config/grade_a_literature_screen.csv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/hr_submission")
    args = parser.parse_args()

    candidates = build_candidates(read_csv(args.atlas), read_csv(args.perturbation), read_csv(args.literature_screen))
    if len(candidates) != 35:
        raise ValueError(f"Expected 35 binding-plus-perturbation candidates, observed {len(candidates)}")
    fields = list(candidates[0])
    write_csv(args.output_dir / "rin_candidate_evidence_table.csv", candidates, fields)

    summary: dict[str, object] = {
        "resource_version": "hr_submission_bundle_v1",
        "candidate_count": len(candidates),
        "grade_counts": dict(sorted(Counter(str(row["candidate_grade"]) for row in candidates).items())),
        "four_cohort_consensus_candidates": sum(bool(row["four_cohort_direction_consensus"]) for row in candidates),
        "protein_direction_concordant_candidates": sum(bool(row["protein_direction_concordant"]) for row in candidates),
        "differential_protein_candidates": sum(bool(row["protein_differential"]) for row in candidates),
        "differential_phosphosite_candidates": sum(bool(row["phosphosite_differential"]) for row in candidates),
        "lower_in_rin1": sum(row["rin1_response"] == "lower_in_rin1" for row in candidates),
        "higher_in_rin1": sum(row["rin1_response"] == "higher_in_rin1" for row in candidates),
        "annotation_theme_counts": dict(sorted(Counter(str(row["annotation_theme"]) for row in candidates).items())),
        "literature_screened_grade_a_candidates": sum(
            row["candidate_grade"] == "A"
            and literature_screen_complete(row["literature_screen_status"])
            for row in candidates
        ),
        "literature_pending_grade_a_candidates": sum(
            row["candidate_grade"] == "A"
            and not literature_screen_complete(row["literature_screen_status"])
            for row in candidates
        ),
        "grading_excludes_studyshield": True,
        "claim_boundary": "Evidence grades prioritise candidates and do not establish direct causality, regulatory sign, or mechanism.",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "rin_candidate_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    protein = json.loads((ROOT / "results/PXD051570_summary.json").read_text())
    binding = json.loads((ROOT / "results/GSE40257_rin_binding_summary.json").read_text())
    perturbation = json.loads((ROOT / "results/GSE210589/deseq2/GSE210589_deseq2_summary.json").read_text())
    tissue = json.loads((ROOT / "results/SRP109982_summary.json").read_text())
    external = read_csv(ROOT / "results/external_validation_GSE78733/GSE78733_frozen_external_benchmark.csv")
    go_rows = read_csv(ROOT / "results/go_enrichment_propagated_matched_nulls.csv")
    figure_dir = args.output_dir / "figures"
    figure1_design(summary, figure_dir / "figure1_study_design.svg")
    figure2_transfer(ROOT, external, figure_dir / "figure2_expression_transfer.svg")
    figure3_orthogonal(protein, binding, perturbation, figure_dir / "figure3_orthogonal_support.svg")
    figure4_applicability(tissue, ROOT, figure_dir / "figure4_tissue_applicability.svg")
    figure5_modules(candidates, go_rows, figure_dir / "figure5_functional_modules.svg")
    figure6_candidates(candidates, figure_dir / "figure6_candidate_landscape.svg")
    print(f"Wrote HR submission bundle with {len(candidates)} candidates and six figures to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
