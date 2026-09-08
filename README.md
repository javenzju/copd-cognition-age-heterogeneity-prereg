# Lung function → cognitive function: two-sample MR — preregistration & reproducible code

This repository contains the preregistration protocol and the reproducible analysis code for:

> *Genetically predicted lower lung function and poorer cognitive function: a two-sample Mendelian randomization study*

## Preregistration

- `preregistration_protocol.md` — original analysis plan (retrospective preregistration).
- `MR_redesign_protocol.md` — revised plan after strategic refocusing of the manuscript on the MR headliner (lung function → cognition) with the individual-level HRS analysis demoted to an exploratory, hypothesis-generating contrast.

## Analysis code (OpenGWAS / TwoSampleMR, R 4.5.3)

- `code/mr_run.R` — main two-sample MR (FEV1, FEV1%pred, COPD liability → intelligence, cognitive performance, Alzheimer's disease, hippocampal volume).
- `code/mr_consolidate.R` — recomputes honest sample sizes, MR-Egger intercept, Cochran Q, manual Steiger directionality, and leave-one-out from cached harmonised data.
- `code/mr_mvmr.R`, `code/mr_mvmr_cog.R` — multivariable MR of FEV1 and smoking initiation (GSCAN `ieu-b-4877`) on cognition.
- `code/mr_bidirectional.R` — reverse-direction MR (intelligence → FEV1%pred).
- `code/mr_overlap_check.R` — sample-overlap diagnostics between exposure and outcome GWAS.

## Results tables

- `results/mr_layer1.csv` — primary Layer-1 MR estimates.
- `results/mr_robust_final.csv` — 15-pair robustness panel (Egger, weighted median/mode, Cochran Q, leave-one-out, manual Steiger).

## Instrumentation note

The correct GSCAN smoking-initiation instrument is `ieu-b-4877`. The legacy `ieu-a-1129` now resolves to breast cancer in OpenGWAS and must **not** be used.

## Data availability

GWAS summary statistics: OpenGWAS (https://gwas.mrcieu.ac.uk/). This repository is archived and citable via Zenodo: **DOI 10.5281/zenodo.22663699** (https://doi.org/10.5281/zenodo.22663699).
