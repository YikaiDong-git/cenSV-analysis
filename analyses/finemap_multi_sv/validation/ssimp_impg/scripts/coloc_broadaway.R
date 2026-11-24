suppressPackageStartupMessages(library(coloc))
P <- "${CENSV_WORKSPACE}"
WD <- file.path(P, "Project1_Centromere/analyses/finemap_multi_sv/validation/ssimp_impg/work/coloc_pros1")
cat(sprintf("coloc pkg version: %s\n", as.character(packageVersion("coloc"))))
NEQTL <- 1183L
run <- function(study, s_cc) {
  fp <- file.path(WD, paste0("coloc_input_broadaway_", study, ".tsv"))
  if (!file.exists(fp)) { cat(sprintf("\n[%s] input missing\n", study)); return(invisible()) }
  df <- read.table(fp, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df <- df[is.finite(df$beta_eqtl) & df$se_eqtl > 0 & is.finite(df$beta_gwas) & df$se_gwas > 0 &
           df$maf_eqtl > 0 & df$maf_eqtl < 1 & df$maf_gwas > 0 & df$maf_gwas < 1 & !duplicated(df$snp), ]
  if (nrow(df) < 5) { cat(sprintf("\n[%s] too few SNPs (%d)\n", study, nrow(df))); return(invisible()) }
  D_eqtl <- list(beta = df$beta_eqtl, varbeta = df$se_eqtl^2, snp = df$snp, type = "quant",
                 N = NEQTL, MAF = df$maf_eqtl, sdY = 1)
  D_gwas <- list(beta = df$beta_gwas, varbeta = df$se_gwas^2, snp = df$snp, type = "cc",
                 N = max(df$n_gwas), MAF = df$maf_gwas, s = s_cc)
  zg <- df$beta_gwas / df$se_gwas; ze <- df$beta_eqtl / df$se_eqtl
  cat(sprintf("\n===== %s x Broadaway-Liver-PROS1 (nsnps=%d) =====\n", study, nrow(df)))
  cat(sprintf("  VTE max|z|=%.2f @ %d ; eQTL max|z|=%.2f @ %d\n",
              max(abs(zg)), df$pos[which.max(abs(zg))], max(abs(ze)), df$pos[which.max(abs(ze))]))
  res <- coloc.abf(D_eqtl, D_gwas); s <- res$summary
  cat(sprintf("  PP.H0=%.3f  PP.H1=%.3f  PP.H2=%.3f  PP.H3=%.3f  PP.H4=%.3f\n",
              s["PP.H0.abf"], s["PP.H1.abf"], s["PP.H2.abf"], s["PP.H3.abf"], s["PP.H4.abf"]))
  cat(sprintf("  H4/(H3+H4)=%.3f\n", s["PP.H4.abf"] / (s["PP.H3.abf"] + s["PP.H4.abf"])))
  o <- res$results[order(-res$results$SNP.PP.H4), c("snp", "SNP.PP.H4")]
  cat("  top 5 SNPs by SNP.PP.H4:\n"); print(head(o, 5), row.names = FALSE)
}
run("finngen", 0.05)
run("ukbb_neale", 0.02)
cat("\nDONE_COLOC_BROADAWAY\n")
