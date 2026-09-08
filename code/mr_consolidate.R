#!/usr/bin/env Rscript
# mr_consolidate.R  (v2 -- fixed Steiger + LOO extraction)
# Authoritative consolidation of MR Layer-1 results from mr_results.rds.
#  - analysis n pulled from each pair's `res` (mr() drops SNPs with non-positive
#    variance / missing) -> honest instrument count for the manuscript.
#  - Steiger directionality via MendelianRandomization::directionality_test
#    (requires eaf + sample sizes, passed from cached har).
#  - Leave-one-out from cached har (no method column in LOO output -> filter SNP!="All").
#  - Egger intercept p + Cochran Q p from stored het/plei.
# Outputs: mr_layer1.csv (corrected main table), mr_robust_final.csv (robustness).

suppressWarnings({ library(TwoSampleMR); library(MendelianRandomization); library(ieugwasr) })
options(ieugwasr_api = "https://api.opengwas.io/api/")
options(timeout = 600)
tok <- Sys.getenv("OPENGWAS_JWT"); if (nzchar(tok)) Sys.setenv(OPENGWAS_JWT = tok)

r <- readRDS("mr_results.rds")
L1 <- r$layer1
pair_names <- names(L1)

old <- read.csv("mr_layer1.csv", stringsAsFactors = FALSE)  # preserve F_stat from prior run
fmap <- setNames(old$F_stat, old$pair)

methods <- c("Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode")
rows_main <- list(); rows_rob <- list()
con <- file("mr_consolidate_summary.txt", "w")
cat("==== MR consolidation v2 ====\n", file = con)

for (pn in pair_names) {
  e <- L1[[pn]]; har <- e$har; res <- e$res
  het <- if ("het" %in% names(e)) e$het else tryCatch(mr_heterogeneity(har), error = function(x) NULL)
  plei <- if ("plei" %in% names(e)) e$plei else tryCatch(mr_pleiotropy_test(har), error = function(x) NULL)
  analysis_n <- if (!is.null(res) && nrow(res) > 0) as.integer(res$nsnp[1]) else nrow(har)
  binary <- e$binary
  Fstat <- if (is.null(fmap[[pn]])) NA else fmap[[pn]]

  for (m in methods) {
    ri <- which(res$method == m)[1]
    if (length(ri) == 0 || is.na(ri)) next
    b <- res$b[ri]; se <- res$se[ri]; p <- res$pval[ri]
    OR <- OR_lo <- OR_hi <- NA
    if (binary) { OR <- exp(b); OR_lo <- exp(b - 1.96 * se); OR_hi <- exp(b + 1.96 * se) }
    rows_main[[length(rows_main) + 1]] <- data.frame(
      analysis = "Layer1", pair = pn, outcome_binary = binary, method = m,
      b = b, se = se, p = p, OR = OR, OR_lo = OR_lo, OR_hi = OR_hi,
      n_snps = analysis_n, F_stat = Fstat,
      egger_intercept_p = if (!is.null(plei)) plei$pval else NA,
      cochran_Q_p = if (!is.null(het)) het$Q_pval[which(het$method == m)[1]] else NA,
      stringsAsFactors = FALSE)
  }

  # ---- Steiger directionality (variance-explained test, manual/transparent) ----
  # Compares the variance in exposure vs outcome explained by the IVs.
  # If R2_exposure > R2_outcome, direction exposure->outcome is supported.
  steiger_p <- steiger_rho <- steiger_dir <- steiger_R2x <- steiger_R2y <- NA
  keep <- har$mr_keep & !is.na(har$beta.exposure) & !is.na(har$beta.outcome) &
         !is.na(har$se.exposure) & !is.na(har$se.outcome) & !is.na(har$eaf.exposure)
  if (sum(keep, na.rm = TRUE) >= 3) {
    hk <- har[keep, ]
    maf <- pmin(hk$eaf.exposure, 1 - hk$eaf.exposure)
    vx <- 2 * maf * (1 - maf) * hk$beta.exposure^2
    vy <- 2 * maf * (1 - maf) * hk$beta.outcome^2
    se_vx <- 4 * maf * (1 - maf) * abs(hk$beta.exposure) * hk$se.exposure
    se_vy <- 4 * maf * (1 - maf) * abs(hk$beta.outcome) * hk$se.outcome
    R2x <- mean(vx, na.rm = TRUE); R2y <- mean(vy, na.rm = TRUE)
    var_R2x <- mean(se_vx^2, na.rm = TRUE) / sum(keep, na.rm = TRUE)
    var_R2y <- mean(se_vy^2, na.rm = TRUE) / sum(keep, na.rm = TRUE)
    z <- (R2x - R2y) / sqrt(var_R2x + var_R2y)
    p <- 2 * pnorm(-abs(z))
    steiger_R2x <- R2x; steiger_R2y <- R2y
    steiger_rho <- z; steiger_p <- p
    steiger_dir <- ifelse(R2x > R2y, "to_outcome", "to_exposure")
    cat(sprintf("  [%s] Steiger R2_exp=%.3g R2_out=%.3g z=%.2f p=%.2g dir=%s\n",
                pn, R2x, R2y, z, p, steiger_dir), file = con)
  }

  # ---- Leave-one-out (only for significant IVW pairs) ----
  loo_min <- loo_max <- loo_all_b <- loo_all_p <- NA
  ivw_p <- if (!is.null(res)) res$pval[which(res$method == "Inverse variance weighted")[1]] else NA
  if (!is.na(ivw_p) && ivw_p < 0.05) {
    loo <- tryCatch(mr_leaveoneout(har), error = function(x) { cat(sprintf("  [%s] LOO err: %s\n", pn, conditionMessage(x)), file = con); NULL })
    if (!is.null(loo)) {
      per <- loo[loo$SNP != "All", ]
      al <- loo[loo$SNP == "All", ]
      if (nrow(per) > 0) { loo_min <- min(per$b, na.rm = TRUE); loo_max <- max(per$b, na.rm = TRUE) }
      if (nrow(al) > 0) { loo_all_b <- al$b; loo_all_p <- al$p }
      cat(sprintf("  [%s] LOO IVW b=[%.4f, %.4f] (All b=%.4f p=%.3g)\n", pn, loo_min, loo_max, loo_all_b, loo_all_p), file = con)
    }
  }

  rows_rob[[length(rows_rob) + 1]] <- data.frame(
    pair = pn, analysis_n = analysis_n,
    egger_intercept_p = if (!is.null(plei)) plei$pval else NA,
    cochran_Q_p = if (!is.null(het)) het$Q_pval[which(het$method == "Inverse variance weighted")[1]] else NA,
    steiger_R2_exp = steiger_R2x, steiger_R2_out = steiger_R2y,
    steiger_z = steiger_rho, steiger_p = steiger_p, steiger_dir = steiger_dir,
    loo_min_b = loo_min, loo_max_b = loo_max, loo_all_b = loo_all_b, loo_all_p = loo_all_p,
    stringsAsFactors = FALSE)
}
close(con)

main_df <- do.call(rbind, rows_main)
rob_df <- do.call(rbind, rows_rob)
write.csv(main_df, "mr_layer1.csv", row.names = FALSE)
write.csv(rob_df, "mr_robust_final.csv", row.names = FALSE)
cat(sprintf("Wrote mr_layer1.csv (%d rows) and mr_robust_final.csv (%d rows)\n", nrow(main_df), nrow(rob_df)))
