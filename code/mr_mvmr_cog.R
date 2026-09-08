# Manual MVMR ({FEV1, smoking} -> CogPerf) using the COMMON instrument set.
# Key fix: mv_extract_exposures(c(FEV1, smoking)) returns the joint instrument set
# (each SNP instruments BOTH exposures); separate extract_instruments give disjoint sets (0 overlap).
# Manual MVMR-IVW estimator (Burgess & Thompson 2015):
#   Y = X1*b1 + X2*b2 + e,  Var(e_j) = seY^2 + (X1*seX1)^2 + (X2*seX2)^2
suppressWarnings({ library(TwoSampleMR); library(ieugwasr) })
options(ieugwasr_api = "https://api.opengwas.io/api/")
options(timeout = 600)

# ---- Step 1: common instruments (crash-prone proxy lookup; cached across runs) ----
if (!file.exists("exps_common.rds")) {
  cat("extracting common instruments (mv_extract_exposures)...\n"); flush(stdout())
  exps <- mv_extract_exposures(c("ukb-b-19657", "ieu-b-4877"), pval_threshold = 5e-8)
  saveRDS(exps, "exps_common.rds")
  cat(sprintf("common instrument SNPs: %d\n", nrow(exps)), file = stderr())
} else {
  exps <- readRDS("exps_common.rds")
  cat(sprintf("loaded %d common SNPs from cache\n", nrow(exps)), file = stderr())
}
snps <- unique(exps$SNP)
cat(sprintf("unique SNPs: %d\n", length(snps)), file = stderr())

# ---- Step 2: CogPerf outcome for those SNPs (cached) ----
if (!file.exists("out_cog2.rds")) {
  out <- extract_outcome_data(snps = snps, outcomes = "ebi-a-GCST006572", proxies = FALSE)
  saveRDS(out, "out_cog2.rds")
} else out <- readRDS("out_cog2.rds")
out$SNP <- as.character(out$SNP)
out$effect_allele.outcome <- as.character(out$effect_allele.outcome)
out$other_allele.outcome  <- as.character(out$other_allele.outcome)

# ---- Step 3: build per-SNP exposure effects from common set ----
fe_sub <- exps[exps$id.exposure == "ukb-b-19657", ]
sm_sub <- exps[exps$id.exposure == "ieu-b-4877",  ]
for (df in c("fe_sub", "sm_sub")) {
  assign(df, within(get(df), { SNP <- as.character(SNP); effect_allele.exposure <- as.character(effect_allele.exposure) }))
}
comp <- c(A = "T", T = "A", C = "G", G = "C")
align_to_outcome <- function(ea, rea, roa, beta) {
  if (ea == rea) return(beta)
  if (ea == roa) return(-beta)
  if (!is.na(comp[ea]) && comp[ea] == rea) return(-beta)
  if (!is.na(comp[ea]) && comp[ea] == roa) return(beta)
  NA
}

rows <- list()
for (s in snps) {
  o <- out[out$SNP == s, ]; if (nrow(o) == 0) next; o <- o[1, ]
  f <- fe_sub[fe_sub$SNP == s, ]; if (nrow(f) == 0) next; f <- f[1, ]
  m <- sm_sub[sm_sub$SNP == s, ]; if (nrow(m) == 0) next; m <- m[1, ]
  bf <- align_to_outcome(f$effect_allele.exposure, o$effect_allele.outcome, o$other_allele.outcome, f$beta.exposure)
  bs <- align_to_outcome(m$effect_allele.exposure, o$effect_allele.outcome, o$other_allele.outcome, m$beta.exposure)
  if (is.na(bf) || is.na(bs)) next
  rows <- append(rows, list(list(Y = o$beta.outcome, seY = o$se.outcome,
                                  X1 = bf, seX1 = f$se.exposure, X2 = bs, seX2 = m$se.exposure)))
}
d <- do.call(rbind, lapply(rows, as.data.frame))
cat(sprintf("aligned SNPs for CogPerf MVMR: %d\n", nrow(d)), file = stderr())

Y <- d$Y; X <- cbind(d$X1, d$X2)
W <- 1 / (d$seY^2 + (d$X1 * d$seX1)^2 + (d$X2 * d$seX2)^2)
XtW <- t(X) %*% (W * X); XtwY <- t(X) %*% (W * Y)
beta <- solve(XtW) %*% XtwY; covb <- solve(XtW)
se <- sqrt(diag(covb)); z <- beta / se; p <- 2 * pnorm(-abs(z))

con <- file("mr_mvmr_cog.txt", "w")
cat("===== MVMR: {FEV1, smoking} -> ebi-a-GCST006572 (CogPerf) =====\n", file = con)
cat(sprintf("aligned SNPs: %d\n", nrow(d)), file = con)
cat(sprintf("  -> FEV1 direct effect (adj smoking):    b=%.4f se=%.4f z=%.3f p=%.3g\n", beta[1], se[1], z[1], p[1]), file = con)
cat(sprintf("  -> Smoking direct effect (adj FEV1):    b=%.4f se=%.4f z=%.3f p=%.3g\n", beta[2], se[2], z[2], p[2]), file = con)
close(con)
cat("CogPerf MVMR written to mr_mvmr_cog.txt\n")
print(data.frame(exposure = c("FEV1(ukb-b-19657)", "smoking(ieu-b-4877)"),
                 b = as.numeric(beta), se = as.numeric(se), z = as.numeric(z), p = as.numeric(p)))
