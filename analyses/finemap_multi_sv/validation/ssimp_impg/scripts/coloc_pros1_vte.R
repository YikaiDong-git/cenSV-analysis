suppressPackageStartupMessages(library(coloc))
P <- "${CENSV_WORKSPACE}"
WD <- file.path(P, "Project1_Centromere/analyses/finemap_multi_sv/validation/ssimp_impg/work/coloc_pros1")
cat(sprintf("coloc pkg version: %s\n", as.character(packageVersion("coloc"))))

run <- function(study, s_cc) {
  fp <- file.path(WD, paste0("coloc_input_", study, ".tsv"))
  if (!file.exists(fp)) { cat(sprintf("\n[%s] input missing\n", study)); return(invisible()) }
  df <- read.table(fp, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df <- df[is.finite(df$z_eqtl) & is.finite(df$beta_gwas) & df$se_gwas > 0 &
           df$maf > 0 & df$maf < 1 & !duplicated(df$snp), ]
  if (nrow(df) < 5) { cat(sprintf("\n[%s] too few SNPs (%d)\n", study, nrow(df))); return(invisible()) }
  v <- 1 / (2 * df$maf * (1 - df$maf) * (df$n_eqtl + df$z_eqtl^2))
  D_eqtl <- list(beta = df$z_eqtl * sqrt(v), varbeta = v, snp = df$snp,
                 type = "quant", N = max(df$n_eqtl), MAF = df$maf, sdY = 1)
  D_gwas <- list(beta = df$beta_gwas, varbeta = df$se_gwas^2, snp = df$snp,
                 type = "cc", N = max(df$n_gwas), MAF = df$maf, s = s_cc)
  cat(sprintf("\n================= %s x eQTLGen-PROS1 (nsnps=%d) =================\n", study, nrow(df)))
  cat(sprintf("  GWAS region: max|z|=%.2f @ pos=%d ; eQTL max|z|=%.2f @ pos=%d\n",
              max(abs(df$beta_gwas / df$se_gwas)), df$pos[which.max(abs(df$beta_gwas / df$se_gwas))],
              max(abs(df$z_eqtl)), df$pos[which.max(abs(df$z_eqtl))]))
  res <- coloc.abf(D_eqtl, D_gwas)
  s <- res$summary
  cat(sprintf("  PP.H0=%.3f  PP.H1=%.3f  PP.H2=%.3f  PP.H3=%.3f  PP.H4=%.3f\n",
              s["PP.H0.abf"], s["PP.H1.abf"], s["PP.H2.abf"], s["PP.H3.abf"], s["PP.H4.abf"]))
  cat(sprintf("  H4/(H3+H4) = %.3f  (conditional prob of shared causal | a signal in both)\n",
              s["PP.H4.abf"] / (s["PP.H3.abf"] + s["PP.H4.abf"])))
  o <- res$results[order(-res$results$SNP.PP.H4), c("snp", "SNP.PP.H4")]
  cat("  top 5 SNPs by SNP.PP.H4:\n"); print(head(o, 5), row.names = FALSE)
}

run("finngen", 0.05)       # FinnGen I9_VTE (powered, Finnish)
run("ukbb_neale", 0.02)    # UKBB Neale I9_VTE (EUR, underpowered; ancestry-matched to eQTLGen)
cat("\nDONE_COLOC\n")
