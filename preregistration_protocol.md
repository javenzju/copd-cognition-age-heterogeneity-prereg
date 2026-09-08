# Preregistration / Transparent Analysis Plan
## Age-Dependent Heterogeneity in the COPD–Cognitive Decline Association

**Status:** This document is published as a public, timestamped analysis plan for the manuscript
*"Age-Dependent Heterogeneity in the Association Between COPD and Cognitive Decline: Longitudinal Evidence From the Health and Retirement Study."*
It is **NOT a preprint** (no manuscript text is hosted here, per the author's explicit instruction).
It records the analysis plan, the pre-specified decisions, and the results as obtained, transparently.

**Repository:** public GitHub repo (analysis code + aggregate result tables + this protocol).
**Cohort:** Harmonized HRS panel, waves 3–13 (2006–2016), `hhidpn` key.
**Primary outcome:** `cog27` (TICS-based 27-point composite). **Sensitivity:** `cogtot`.

---

## 1. Background & Rationale

Vivek et al. (medRxiv 2025.05.29.25328529; Sci Rep 2026;16:17461, PMID 41981149) established that
*impaired lung function → elevated AD blood biomarkers (NfL, pTau181) → incident dementia* in HRS.
That paper occupies the "lung function → AD-pathology → dementia" causal chain.

**This study asks a different question:** *Does the magnitude/direction of the COPD–cognitive-decline
association itself vary by age at onset, and can selective survival or biologic aging explain that variation?*
This is a heterogeneity/biased-estimation question, not a lung→biomarker mediation question, and is not preempted.

## 2. Research Question

Is the association between COPD and cognitive decline **age-dependent**, and is that age-dependence
explainable by (a) selective survival (competing mortality) or (b) biologic-aging markers (DunedinPACE, telomere)?

## 3. Hypotheses (direction-open)

- **H1 (heterogeneity).** The COPD–cognition association differs by baseline age group (<59 vs ≥59).
  Direction not pre-specified.
- **H2 (survival).** Any older-onset attenuation is *not* fully explained by selective survival.
  Tested via IPCW + EMM bound (require excess-mortality multiplier EMM to null the effect).
- **H3 (biologic aging, exploratory, direction-open).** DunedinPACE (but not telomere) partially explains
  the heterogeneity. Tested as COPD×age×pace three-way interaction.
- **H4 (mediation, exploratory, direction-open).** The NfL/pTau181/Hcy pathway differs by age group
  (joint mediation estimated separately in each group; no direction pre-specified).

## 4. Design & Data

- Longitudinal, individual-level two-way fixed-effects (TWFE) regression; HRS 2006–2016 (waves 3–13).
- Exposure: COPD (self-reported ever-COPD + wave-13 physician-diagnosed `copd_wave`).
- Outcome: `cog27` (primary), `cogtot` (sensitivity).
- **Reconciled covariate specification (FULL, used for ALL primary models):**
  `agey_b` + female + raedyrs + race_white + race_black + hispanic + married + smoken + cesd
  + bmi + hibpe + diabe + hearte + stroke + cancre.
  *(This single spec resolves the prior +0.344 vs +0.107 older-group discrepancy between the MIN/09C
  model and the FULL/Stage-14 model — the 3× difference was purely covariate-set artifact.)*

## 5. Effect Modifiers & Split

- **Age group:** primary split at baseline age **59** (matches manuscript history); sensitivity at
  empirical median **61**.
- **DunedinPACE:** `paceofaging` (from HRSPoA_Shared.dta, key `hhidpn`).
- **Telomere:** `telomere_short` / `telomere_log` (already in panel). **Three-way telomere interaction
  DROPPED a priori for power** (COPD∩telomere only ~22–36 persons per age cell; would yield uninformative CIs).
  Telomere retained only as a full-sample (non-stratified) effect modifier.

## 6. Statistical Analysis (pre-specified)

- TWFE with alternating-projection demeaning; SE clustered by `hhidpn`.
- Between-group difference: **Wald test** on (β_younger − β_older).
- Selective survival: IPCW weights + EMM bound (effect nulled only if EMM ≥ observed excess-mortality HR).
- Biologic aging: COPD×age×pace three-way interaction (age<59 and median-split); report MDβ80 (power 80%).
- Mediation: OLS path a/b, NfL/pTau181/Hcy, joint index; bootstrap 2000; run separately by age group.

## 7. Pre-specified Decision / Stopping Rules

- **Title wording:** "age-dependent heterogeneity" if FULL-spec between-group Wald p ≥ 0.05;
  "reversal" only if Wald p < 0.05 in FULL spec. (Decision: **heterogeneity** — Wald p = 0.187 at 59,
  0.052 at empirical median; not a clean reversal.)
- **Three-way null interpretation:** if MDβ80 ≈ main-effect magnitude (underpowered), report as
  *underpowered / uninformative*, NOT as evidence of no moderation. (MDβ80 ≈ 0.53 → underpowered.)
- **Mediation:** report both groups direction-open; do NOT claim "younger stronger" — the data show the
  opposite (older-group pathway significant).

## 8. Results As Obtained (transparency)

| Analysis | Key result |
|---|---|
| Age split (cog27, FULL, 59) | younger β=−0.144 (p=.111); older β=+0.015 (p=.850); Wald p=.187 |
| Age split (cog27, FULL, median 61) | younger β=−0.164 (p=.045); older β=+0.070 (p=.428); Wald p=.052 |
| Selective survival (Stage 9C) | IPCW mean 0.927; EMM≥4 needed; real HR 1.3–2.0 → insufficient |
| COPD×age×pace three-way | null (p=.79–.93); MDβ80≈0.53 (underpowered) |
| Telomere full-sample mod. | copd×telomere_short β=−0.527 (p=.058, borderline) |
| Joint mediation older (≥59) | JOINT=−0.131 (p=.018); NfL p=.017, pTau181 p=.046 |
| Joint mediation younger (<59) | JOINT=−0.033 (p=.50); all n.s. |

## 9. Registration Note

This is a transparent/retrospective analysis plan (analyses were conducted before formal registration,
consistent with the author's decision not to post a preprint). It is provided for methodological
transparency and to lock the analysis decisions publicly. OSF-proper registration would require a
separate OSF login and can be added later if desired.

**Archived snapshot (citable):** this protocol and its companion analysis code are archived with a
Zenodo DOI: **10.5281/zenodo.22663699** (https://doi.org/10.5281/zenodo.22663699).
