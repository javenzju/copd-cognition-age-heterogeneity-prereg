# Sample-overlap check for the MR causal layer.
# Exposure ukb-b-19657 (FEV1 % predicted) is UK-Biobank-only.
# Cognition outcomes (Sniekers/Savage 2018) are UKB-inclusive meta-GWAS -> real overlap.
suppressWarnings({ library(TwoSampleMR); library(ieugwasr) })
options(ieugwasr_api = "https://api.opengwas.io/api/")
options(timeout = 600)

con <- file("mr_overlap_check.txt", "w")
cat("===== SAMPLE OVERLAP CHECK (Exposure = UK Biobank-only) =====\n", file = con)

exp_ids <- c(
  "ukb-b-19657" = "FEV1 %predicted (UK Biobank only)",
  "ukb-b-5113"  = "FEV1 adjusted-for-sex (UK Biobank only)"
)
out_ids <- c(
  "ebi-a-GCST006250" = "Intelligence (Sniekers/Savage 2018)",
  "ebi-a-GCST006572" = "Cognitive performance (Savage 2018)",
  "ebi-a-GCST90027158" = "AD overall (Bellenguez 2022)",
  "ebi-a-GCST002245" = "AD late-onset (Lambert 2013)",
  "ubm-b-2868" = "Right hippocampal volume (UK Biobank)"
)

cat("\n--- Exposure GWAS info ---\n", file = con)
for (id in names(exp_ids)) {
  gi <- tryCatch(gwasinfo(id), error = function(e) NULL)
  if (!is.null(gi)) cat(sprintf("[EXP] %s | %s | n=%s | build=%s\n",
                                id, gi$trait, gi$ncase + gi$ncontrol, gi$build), file = con)
}

cat("\n--- Outcome GWAS info ---\n", file = con)
for (id in names(out_ids)) {
  gi <- tryCatch(gwasinfo(id), error = function(e) NULL)
  if (!is.null(gi)) cat(sprintf("[OUT] %s | %s | n=%s | build=%s\n",
                                id, gi$trait, gi$ncase + gi$ncontrol, gi$build), file = con)
}

# Qualitative overlap classification based on consortium sample composition.
notes <- c(
  "ebi-a-GCST006250" = "UKB-INCLUSIVE meta-GWAS (Sniekers 2018 combines UKB + consortium). Exposure ukb-b-19657 is pure UKB -> SUBSTANTIAL SAMPLE OVERLAP.",
  "ebi-a-GCST006572" = "UKB-INCLUSIVE (Savage 2018). -> SUBSTANTIAL SAMPLE OVERLAP.",
  "ebi-a-GCST90027158" = "IGAP/GR@P case-control, minimal UKB -> LOW overlap.",
  "ebi-a-GCST002245" = "IGAP, minimal UKB -> LOW overlap.",
  "ubm-b-2868" = "Same UKB cohort as exposure -> MAXIMAL overlap (effectively single-sample; null result reported)."
)
cat("\n--- OVERLAP CLASSIFICATION ---\n", file = con)
for (id in names(out_ids)) cat(sprintf("%s: %s\n", id, notes[id]), file = con)

# Candidate non-UKB cognition GWAS for sensitivity (lightweight search).
cat("\n--- Candidate non-UKB cognition GWAS (sensitivity option) ---\n", file = con)
ao <- tryCatch(available_outcomes("cognitive performance"), error = function(e) NULL)
if (!is.null(ao) && nrow(ao) > 0) {
  # prefer entries whose consortium/trait does NOT obviously include UKB
  for (i in seq_len(min(nrow(ao), 25))) {
    cat(sprintf("  %s | %s | n=%s | %s\n",
                ao$id[i], ao$trait[i], ao$ncase[i] + ao$ncontrol[i],
                ifelse(is.null(ao$consortium[i]), "", as.character(ao$consortium[i]))), file = con)
  }
} else {
  cat("  (search unavailable)\n", file = con)
}

cat("\n--- MITIGATION STRATEGY (to report) ---\n", file = con)
cat("1. Primary cognition outcomes (Intelligence, CogPerf) have UKB overlap with the FEV1 exposure.\n", file = con)
cat("2. Overlap bias direction: IVW/Egger drift toward the observational association.\n", file = con)
cat("3. Already shown robust: MR-Egger intercept n.s. (no horizontal pleiotropy); weighted median consistent; manual Steiger supports exposure->outcome.\n", file = con)
cat("4. Sensitivity recommended: re-run primary pairs with a non-UKB cognition GWAS if one is suitable (see candidates above).\n", file = con)
close(con)
cat("Overlap check written to mr_overlap_check.txt\n")
