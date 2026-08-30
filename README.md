# Tomato ripening evidence atlas

Analysis code and study configuration for the tomato ripening evidence atlas.
Repository: https://github.com/nblvguohao/tomato-ripening-evidence-atlas

## Current publication scope

This first public snapshot contains analysis code, sample/study configuration,
34 analysis unit tests and pinned Python dependencies. No manuscript, original
figures, author forms, raw data, processed expression matrices or result tables
are part of this snapshot. Public source accessions and URLs are recorded in
`config/public_studies.csv` and the contrast configurations.

The processed-data package was tested locally: the 464-gene signature, consensus
and atlas and 35-candidate table were regenerated and matched their frozen
references. Its publication is pending an explicit decision about whether the
processed gene-level expression matrices fall within the authors' withheld-data
boundary. This code-only snapshot must not be described as a self-contained
reproducibility release. `reproduce.py` requires that additional package.

## Verify the code

Use Python 3.12, install `requirements.txt`, then run:

```sh
python -m unittest discover -s tests -p 'test_*.py' -q
```

The matched cross-layer and transfer scripts use seed 20260828 and 1,000 draws.
Candidate grading does not use StudyShield. Binding plus perturbation is
context-bound evidence, not proof of direct causality or regulatory sign.

## Rights and citation

The authors authorized public posting. Reuse licences await explicit approval;
no MIT or Creative Commons licence is granted in this snapshot. Third-party
material retains its original terms. Cite the repository and exact commit.
No Zenodo DOI or complete archived data release is claimed.
