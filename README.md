# COPD × Cognitive Decline — Age-Dependent Heterogeneity (Public Analysis Plan / Preregistration)

This repository is a **public, timestamped analysis plan and reproducible-code record** for the manuscript:

> *Age-Dependent Heterogeneity in the Association Between COPD and Cognitive Decline: Longitudinal Evidence From the Health and Retirement Study*

**This is NOT a preprint.** No manuscript text is hosted here (the author chose not to post a preprint).
It records the pre-specified analysis decisions, the analysis code, and the aggregate result tables so the
methodology is transparent and the analysis choices are locked publicly.

## Why this exists
Vivek et al. (Sci Rep 2026;16:17461, medRxiv 2025.05.29.25328529) occupy the
"impaired lung function → AD blood biomarkers → incident dementia" causal chain. This study asks a
**different** question: whether the COPD–cognitive-decline association is *age-dependent*, and whether
selective survival or biologic aging can explain that variation. The angle is distinct and was confirmed
not preempted by a medRxiv/bioRxiv/AAIC search (2026-09-07).

## Contents
- `preregistration_protocol.md` — pre-specified hypotheses, design, covariate-spec reconciliation, decision rules, results-as-obtained.
- `manuscript_JAGS_age_heterogeneity.docx` / `.md` — the submission-ready JAGS (Clinical Investigation) draft.
- `*.py` — reproducible analysis code (run with `E:/program/python.exe` or any pandas/numpy/statsmodels env).
- `*_valid.csv`, `stage15A_threeway.csv`, `stage15B_age_mediation.csv`, `telomere_full_modification.csv` — aggregate results.
- `最终整合报告_年龄异质性方案.md` — integrated report (Chinese).

## Key decisions (locked)
- **Title wording:** "age-dependent heterogeneity" (between-group Wald p = .187 at 59 / .052 at empirical median 61), NOT "reversal".
- **Covariate spec:** single FULL specification (with `agey_b`) for all primary models — resolves the prior +0.344 vs +0.107 older-group discrepancy.
- **Telomere three-way interaction dropped for power**; telomere retained only as full-sample modifier.
- **DunedinPACE three-way null = underpowered** (MDβ80 ≈ 0.53), reported as uninformative, not as "no moderation".
- **Mediation reported direction-open**; the measured NfL/pTau181 pathway mediates only the older-onset effect.

## Data
Source: Harmonized HRS (waves 3–13, 2006–2016), `hhidpn` key. Individual-level data are not included; only
aggregate result tables and analysis code are shared, consistent with HRS data-use agreements.
