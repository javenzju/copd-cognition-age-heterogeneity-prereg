# Multivariable MR: FEV1 direct effect on cognition, controlling for smoking initiation.
# Goal: defend against "FEV1->cognition is just smoking pleiotropy" by conditioning on a smoking instrument.
# IMPORTANT: ieu-a-1129 is BREAST CANCER (BCAC) in the current OpenGWAS rebuild -> NOT smoking.
#            Correct GSCAN smoking-initiation ID = ieu-b-4877 (verified: GSCAN, n=607291, 2019).
suppressWarnings({ library(TwoSampleMR); library(ieugwasr) })
options(ieugwasr_api = "https://api.opengwas.io/api/")
options(timeout = 600)

con <- file("mr_mvmr.txt", "w")
cat("===== MULTIVARIABLE MR: FEV1 + smoking initiation -> cognition =====\n", file = con)
cat("Exposure 1: ukb-b-19657 (FEV1 % predicted; UK-Biobank-only)\n", file = con)
cat("Exposure 2: ieu-b-4877 (GSCAN smoking initiation; ncase=311629, ncontrol=321173, 2019)\n", file = con)
cat("Note: ieu-a-1129 resolves to Breast cancer (BCAC) in current OpenGWAS; corrected to ieu-b-4877.\n\n", file = con)

exp_ids <- c("ukb-b-19657", "ieu-b-4877")
for (out_id in c("ebi-a-GCST006250", "ebi-a-GCST006572")) {
  cat(sprintf("===== MVMR: {FEV1, smoking} -> %s =====\n", out_id), file = con)

  exps <- tryCatch(
    mv_extract_exposures(exp_ids, pval_threshold = 5e-8),
    error = function(e) { cat("extract err:", conditionMessage(e), "\n", file = con); NULL })
  if (is.null(exps)) { cat("  (no exposures extracted)\n", file = con); next }
  if (!is.data.frame(exps)) exps <- do.call(rbind, exps)   # safety: list -> df
  snps <- unique(exps$SNP)
  cat(sprintf("  extracted %d instrument rows across %d exposures; %d unique SNPs\n",
              nrow(exps), length(unique(exps$id.exposure)), length(snps)), file = con)

  out <- tryCatch(
    extract_outcome_data(snps = snps, outcomes = out_id, proxies = TRUE),
    error = function(e) { cat("outcome err:", conditionMessage(e), "\n", file = con); NULL })
  if (is.null(out)) next

  harm <- tryCatch(
    mv_harmonise_data(exps, out),
    error = function(e) { cat("harmonise err:", conditionMessage(e), "\n", file = con); NULL })
  if (is.null(harm)) next
  cat(sprintf("  harmonised %d SNPs for MVMR\n", nrow(harm)), file = con)

  res <- tryCatch(
    mv_multiple(harm),
    error = function(e) { cat("mvmr err:", conditionMessage(e), "\n", file = con); NULL })
  if (!is.null(res)) {
    # mv_multiple returns a list with a $result data.frame (and possibly other slots)
    res_df <- tryCatch({
      r <- if (is.list(res) && !is.data.frame(res)) res$result else res
      if (is.null(r) || nrow(r) == 0) stop("empty MVMR result")
      r
    }, error = function(e) { cat("  parse err:", conditionMessage(e), "\n", file = con); NULL })
    if (!is.null(res_df)) {
      print(res_df)
      capture.output(print(res_df), file = con)
      fev1 <- tryCatch(res_df[grepl("ukb-b-19657", res_df$exposure), ], error = function(e) NULL)
      smk  <- tryCatch(res_df[grepl("ieu-b-4877", res_df$exposure), ], error = function(e) NULL)
      if (!is.null(fev1) && !is.null(smk)) {
        cat(sprintf("  -> FEV1 direct effect (adjusted for smoking): b=%.4f se=%.4f p=%.3g\n",
                    fev1$b, fev1$se, fev1$pval), file = con)
        cat(sprintf("  -> Smoking direct effect (adjusted for FEV1): b=%.4f se=%.4f p=%.3g\n",
                    smk$b, smk$se, smk$pval), file = con)
      } else cat("  (could not extract per-exposure rows)\n", file = con)
    }
  }
}
close(con)
cat("MVMR written to mr_mvmr.txt\n")
