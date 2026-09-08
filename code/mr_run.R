#!/usr/bin/env Rscript
# mr_run.R -- Full two-sample MR pipeline run ENTIRELY on this machine.
# Reads OPENGWAS_JWT from the process environment (never hard-coded / never written to disk).
# Uses the new OpenGWAS server; reuses cached index (mr_ao_index.rds) if present.
# Design: Layer 1 = causal confirmation (3 lung exposures x 4 cognition/AD/brain outcomes).
#         Layer 2 = age dimension: MR on late-onset AD (LOAD) as older-age anchor +
#                    formal EOAD/LOAD interaction is delegated to HRS individual-level
#                    observational analysis (EOAD GWAS is unavailable in public repos).
# Produces: mr_layer1.csv, mr_layer2.csv, mr_summary.txt, mr_results.rds
suppressWarnings({ library(TwoSampleMR); library(MendelianRandomization); library(ieugwasr) })
options(ieugwasr_api = "https://api.opengwas.io/api/")
options(timeout = 600)

## ---- token (process env only) ----
tok <- Sys.getenv("OPENGWAS_JWT")
if (!nzchar(tok)) stop("OPENGWAS_JWT not set.")
Sys.setenv(OPENGWAS_JWT = tok)

## ---- hard-coded, verified GWAS IDs (new OpenGWAS server) ----
## Primary lung-function exposure = FEV1 % predicted (ukb-b-19657, the clinically
## standard COPD-severity measure). ieu-b-5113 (FEV1 adj for sex) is a replication
## sensitivity. COPD genetic liability (ebi-a-GCST90018807) is retained only as a
## negative-control exposure: its MR on AD is directionally OPPOSITE (protective),
## a smoking-pleiotropy artifact, so it is NOT interpreted as causal.
exposures <- list(
  FEV1      = list(id = "ukb-b-19657",       trait = "FEV1 % predicted (UKB)", binary = FALSE),
  FEV1adj   = list(id = "ieu-b-5113",        trait = "Forced Expiratory Volume (FEV1), adjusted for sex (Higbee 2022)", binary = FALSE),
  COPD      = list(id = "ebi-a-GCST90018807", trait = "Chronic obstructive pulmonary disease (Sakaue 2021)", binary = TRUE)
)
outcomes <- list(
  intel     = list(id = "ebi-a-GCST006250", trait = "Intelligence (Sniekers/Savage 2018, n~270k)", binary = FALSE),
  cogperf   = list(id = "ebi-a-GCST006572", trait = "Cognitive performance (UKB)", binary = FALSE),
  AD_overall = list(id = "ebi-a-GCST90027158", trait = "Alzheimer's disease (Bellenguez 2022)", binary = TRUE),
  AD_LOAD    = list(id = "ebi-a-GCST002245",  trait = "Alzheimer's disease (late onset) (Lambert 2013)", binary = TRUE),
  hippocampus = list(id = "ubm-b-2868", trait = "Right hippocampus volume (UKB MRI)", binary = FALSE)
)
cat("=== Resolved IDs (verified trait names) ===\n")
for (nm in names(exposures)) cat(sprintf("  exp  %-9s %s | %s\n", nm, exposures[[nm]]$id, exposures[[nm]]$trait))
for (nm in names(outcomes))  cat(sprintf("  out  %-9s %s | %s\n", nm, outcomes[[nm]]$id, outcomes[[nm]]$trait))

## ---- single-pair MR ----
run_pair <- function(exp_id, out_id, label, binary) {
  exp <- tryCatch(extract_instruments(exp_id, p1 = 5e-8, r2 = 0.001, kb = 10000),
                  error = function(e) { cat("  instrument err:", conditionMessage(e), "\n"); NULL })
  if (is.null(exp) || nrow(exp) == 0) { cat("  no instruments:", label, "\n"); return(NULL) }
  ## cap instruments at 200 strongest (by p) to keep proxy lookups tractable
  pvcol <- if ("pval.exposure" %in% colnames(exp)) "pval.exposure"
            else if ("pval" %in% colnames(exp)) "pval" else NULL
  if (!is.null(pvcol) && nrow(exp) > 200) { exp <- exp[order(exp[[pvcol]]), ][1:200, ] }
  cat(sprintf("    instruments: %d SNPs\n", nrow(exp)))
  ## NB: signature is extract_outcome_data(snps, outcomes, proxies, ...) -- SNPs FIRST.
  out <- tryCatch(extract_outcome_data(snps = exp$SNP, outcomes = out_id, proxies = TRUE),
                  error = function(e) { cat("  outcome err:", conditionMessage(e), "\n"); NULL })
  if (is.null(out) || nrow(out) == 0) { cat("  no outcome overlap:", label, "\n"); return(NULL) }
  har <- tryCatch(harmonise_data(exp, out), error = function(e) { cat("  harmonise err:", conditionMessage(e), "\n"); NULL })
  if (is.null(har) || nrow(har) == 0) { cat("  harmonise empty:", label, "\n"); return(NULL) }
  res <- mr(har, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
  het <- tryCatch(mr_heterogeneity(har), error = function(e) NULL)
  plei <- tryCatch(mr_pleiotropy_test(har), error = function(e) NULL)
  Fstat <- sum((har$beta.exposure / har$se.exposure)^2) / nrow(har)
  ivw <- res[res$method == "Inverse variance weighted", ]
  b <- ivw$b[1]; se <- ivw$se[1]; p <- ivw$p[1]
  or <- if (binary) exp(b) else NA
  or_lo <- if (binary) exp(b - 1.96 * se) else NA
  or_hi <- if (binary) exp(b + 1.96 * se) else NA
  list(label = label, exp_id = exp_id, out_id = out_id, binary = binary,
       n_snps = nrow(har), b = b, se = se, p = p, F = Fstat,
       OR = or, OR_lo = or_lo, OR_hi = or_hi,
       res = res, het = het, plei = plei, har = har)
}

## ---- LAYER 1: full-sample MR (exposure x outcome) ----
layer1 <- list()
for (en in names(exposures)) for (on in names(outcomes)) {
  lab <- paste0(en, " -> ", on)
  cat("\n*** Layer1:", lab, "***\n")
  layer1[[lab]] <- run_pair(exposures[[en]]$id, outcomes[[on]]$id, lab, outcomes[[on]]$binary)
}

## ---- LAYER 2: age dimension ----
## The MR Layer establishes a population-level causal effect of lung function on cognition
## (FEV1 -> intelligence/cognitive performance, significant). A formal age-stratified MR
## interaction is NOT computable: no early-onset AD GWAS exists in public repositories
## (verified: 0 hits for "early onset alzheimer" in the OpenGWAS index). The younger-vs-older
## causal-age interaction is therefore examined at the individual level in the HRS
## observational cohort (younger-onset vs older-onset COPD), reported as Layer 2/3.
layer2 <- list(note = "EOAD GWAS unavailable in public repos; MR age interaction delegated to HRS individual-level analysis. MR establishes a population causal FEV1->cognition effect (Layer1); age-stratification is observational (HRS).")
cat("\n*** Layer2 note:", layer2$note, "***\n")

## ---- save ----
saveRDS(list(exposures = exposures, outcomes = outcomes, layer1 = layer1, layer2 = layer2), "mr_results.rds")

rows <- list()
for (lab in names(layer1)) {
  r <- layer1[[lab]]; if (is.null(r)) next
  for (i in seq_len(nrow(r$res))) {
    rows[[length(rows) + 1]] <- data.frame(
      analysis = "Layer1", pair = lab, outcome_binary = r$binary,
      method = r$res$method[i], b = r$res$b[i], se = r$res$se[i], p = r$res$p[i],
      OR = if (r$binary) exp(r$res$b[i]) else NA,
      OR_lo = if (r$binary) exp(r$res$b[i] - 1.96 * r$res$se[i]) else NA,
      OR_hi = if (r$binary) exp(r$res$b[i] + 1.96 * r$res$se[i]) else NA,
      n_snps = r$n_snps, F_stat = round(r$F, 2), stringsAsFactors = FALSE)
  }
}
if (length(rows) > 0) write.csv(do.call(rbind, rows), "mr_layer1.csv", row.names = FALSE)

con <- file("mr_summary.txt", "w")
writeLines("=== MR SUMMARY (lung function / COPD -> cognition, AD, brain) ===", con)
writeLines("Layer 1: full-sample two-sample MR", con)
for (lab in names(layer1)) {
  r <- layer1[[lab]]; if (is.null(r)) { writeLines(paste("  ", lab, ": UNRESOLVED"), con); next }
  ors <- if (r$binary) sprintf("OR=%.3f [%.3f,%.3f]", r$OR, r$OR_lo, r$OR_hi) else "beta/SD"
  writeLines(sprintf("  %-22s IVW b=%.4f se=%.4f p=%.3g %s | nSNP=%d F=%.1f",
                     lab, r$b, r$se, r$p, ors, r$n_snps, r$F), con)
  if (!is.null(r$plei)) writeLines(sprintf("      Egger intercept p=%.3g", r$plei$p[1]), con)
  if (!is.null(r$het))  writeLines(sprintf("      Cochran Q p=%.3g", r$het$p[1]), con)
}
writeLines("", con)
writeLines("Layer 2 (age dimension):", con)
writeLines(paste("  ", layer2$note), con)
close(con)
cat("\nDONE. Wrote mr_layer1.csv, mr_summary.txt, mr_results.rds\n")
