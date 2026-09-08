suppressWarnings({library(TwoSampleMR); library(ieugwasr)})
options(ieugwasr_api="https://api.opengwas.io/api/")
con <- file("mr_bidirectional.txt", "w")
cat("=== Bidirectional MR: cognition (intelligence, ebi-a-GCST006250) as exposure -> FEV1percentpred (ukb-b-19657) as outcome ===\n", file=con)

exp <- tryCatch(extract_instruments("ebi-a-GCST006250", p1=5e-8, r2=0.001, kb=10000),
                error=function(e){cat("exp err:",conditionMessage(e),"\n",file=con); NULL})
if(is.null(exp)){ close(con); quit(status=1) }
cat(sprintf("instruments: %d SNPs\n", nrow(exp)), file=con)

out <- tryCatch(extract_outcome_data(snps=exp$SNP, outcomes="ukb-b-19657", proxies=TRUE),
                error=function(e){cat("out err:",conditionMessage(e),"\n",file=con); NULL})
if(is.null(out)){ close(con); quit(status=1) }

har <- tryCatch(harmonise_data(exp, out),
                error=function(e){cat("har err:",conditionMessage(e),"\n",file=con); NULL})
if(is.null(har) || nrow(har)==0){ close(con); quit(status=1) }
cat(sprintf("harmonised: %d SNPs\n", nrow(har)), file=con)

res <- mr(har)
print(res); capture.output(print(res), file=con)
fst <- mean((har$beta.exposure/har$se.exposure)^2, na.rm=TRUE)
cat(sprintf("F (mean of (beta/se)^2 across instruments): %.1f\n", fst), file=con)
cat(sprintf("IVW beta=%.4f se=%.4f p=%.3g\n", res$b[res$method=="Inverse variance weighted"],
            res$se[res$method=="Inverse variance weighted"], res$pval[res$method=="Inverse variance weighted"]), file=con)
cat(sprintf("Egger beta=%.4f se=%.4f p=%.3g\n", res$b[res$method=="MR Egger"],
            res$se[res$method=="MR Egger"], res$pval[res$method=="MR Egger"]), file=con)
cat(sprintf("Weighted median beta=%.4f se=%.4f p=%.3g\n", res$b[res$method=="Weighted median"],
            res$se[res$method=="Weighted median"], res$pval[res$method=="Weighted median"]), file=con)
close(con)
