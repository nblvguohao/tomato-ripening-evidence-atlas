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
`CITATION.cff`. A Zenodo DOI will be recorded only after a real archive has been
published; no DOI is currently claimed.
