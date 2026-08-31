# Tomato ripening evidence atlas — version 1.1.0

Release date: 31 August 2026. Version DOI: https://doi.org/10.5281/zenodo.22185377. This version includes the supplementary analyses and prior-list comparison. Run `python reproduce.py`, `python reproduce_revision.py`, and `python reproduce_prior_comparison.py` after installing requirements.txt. They verify the original four tables, twelve earlier revision tables and four prior-comparison tables.

The prior comparison uses minimal factual gene identifiers, numeric observations and source-row coordinates extracted from Fujisawa 2013 supplements, plus 93 qualifying candidate peak rows. Original article PDF, supplement ZIP/XLS files, descriptions, manuscripts, figures and large raw datasets are excluded. Third-party rights are retained; the project licence does not relicense the original article or supplements. Acquisition and extraction hashes are in data/processed/fujisawa2013/source_provenance.json. Raw-supplement extraction needs the original cited files and xlrd 2.0.2; the portable comparison needs no xlrd.

IMPORTANT occupancy definition: frozen binding/occupancy fields encode gene-associated peaks in at least two independent tables, not necessarily overlapping genomic sites or the original gene-region selection. Solyc04g080540 has no qualifying cross-replicate interval overlap. Candidate grades and membership remain frozen; absent prior-list membership is not novelty. New augmented fields are in prior_expected/candidate_evidence_with_prior.csv; revision_expected retains the previous schema for reproducibility.

## Historical workflow notes (superseded where current definitions differ)

Current machine-readable definitions: data/current_data_dictionary.csv, matching the corrected submission workbook. Legacy revision_expected dictionaries remain byte-preserved for the twelve-table reconstruction and are not the current terminology authority. Coordinate metadata are in the new peak-coordinate table.

### Supplementary scientific analyses

These supplementary analyses are included in release 1.1.0. The original release 1.0.0 remains available at https://doi.org/10.5281/zenodo.22177891.

Install requirements.txt. Run `python reproduce.py` for the original four-table checks and `python reproduce_revision.py` for the revision checks. The latter regenerates five fixed comparisons, exact endpoint-column resampling, candidate contexts, source-selection support, resource-field updates and a post-review worked query. Processed inputs and frozen upstream evidence are supplied. No manuscript, final figures, raw reads, private author metadata or large raw platform file is included. The compact GPL15968 gene universe is traced in revision_input_provenance.json. The new analysis is conditional on normalized inputs and frozen gene sets; it is not a full raw-data reanalysis.

Original phosphosite overlap outputs use a total-protein background and remain historical compatibility outputs only; they are not valid phosphosite-enrichment inference. Revised fields distinguish positive records from absent records with unknown testability. The worked query is an illustration using frozen annotation and existing observations, not a preregistered or validated experimental ranking.

## Original processed-input documentation

# Tomato ripening evidence atlas

Analysis code and compact processed inputs for an evidence-tiered, multi-cohort
tomato ripening atlas. Repository:
https://github.com/nblvguohao/tomato-ripening-evidence-atlas

## Scope and privacy

This repository contains neither the manuscript nor its supplementary prose,
cover letter, author forms, original figure files, private correspondence,
large raw sequencing files, or unpublished follow-on analyses. It starts with
a fresh publication history rather than the working project's history.

Processed gene-level expression matrices, sample labels, frozen evidence tables,
source accessions, analysis scripts and regression tests are included. Retaining
the background genes is necessary for valid controls; providing only the
selected genes would not reproduce those comparisons.

## Reproduce from processed inputs

Use Python 3.12 and the pinned requirements. No API key is needed.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce.py
```

The command verifies file hashes, runs the included analysis tests, reselects
the frozen signature from the two source matrices, reconstructs expression
consensus from independent matrices, rebuilds the evidence graph and atlas,
and reproduces the candidate table. It compares regenerated tables against
the frozen reference tables with explicit numeric tolerances. Results are
written to `reproduced/`; reference inputs are not overwritten.

This is processed-input reproduction, not a claim that all raw-data analyses
can run without downloads. The protein and perturbation results and historical
StudyShield scores are frozen intermediate inputs to the atlas reconstruction.
The original processing scripts are retained for transparency. Additional public
inputs are needed to recompute these upstream layers or their matched nulls.

## Original sources and full upstream recomputation

See `config/public_studies.csv`, `config/cohort_contrasts.csv`, and
`results/source_sha256_manifest.json` for accessions, roles, source URLs and
checksums. Raw data remain at GEO/SRA/ProteomeXchange and the original publishers;
large third-party files are deliberately not redistributed here. Retrieve the
specified versions, verify checksums, and consult each script's `--help` before
raw-data processing. R/DESeq2 is additionally required for the GSE210589 upstream
analysis. The pinned gene annotation is Ensembl Plants release 63; do not silently
substitute a newer annotation or ontology.

The matched cross-layer and transfer scripts use seed **20260828** and 1,000
draws. These seeds must not be replaced with the general pipeline seed.
Candidate grading does not use StudyShield. Binding plus perturbation is
context-bound evidence, not proof of direct causality or regulatory sign.

## Licensing and citation

Public posting of the compact processed matrices and necessary result tables
was explicitly authorized by the authors on 30 August 2026. Original analysis
code is licensed under the MIT License (see `LICENSE`). Author-generated result
tables are licensed under CC BY 4.0 (see `DATA_LICENSE.md`). Third-party data,
annotations and derived portions retain their original terms; these grants do
not relicense third-party rights.

Cite this repository and the exact commit used; author metadata are provided in
`CITATION.cff`. Cite version 1.1.0: https://doi.org/10.5281/zenodo.22185377. The original 1.0.0 archive is retained for version history.
