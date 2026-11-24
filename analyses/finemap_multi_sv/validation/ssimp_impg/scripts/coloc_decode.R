suppressPackageStartupMessages({ library(coloc) })
P <- "${CENSV_WORKSPACE}"
WD <- file.path(P, "Project1_Centromere/analyses/finemap_multi_sv/validation/ssimp_impg/work/coloc_pros1")
cat(sprintf("coloc pkg version: %s\n", as.character(packageVersion("coloc"))))
npf <- file.path(WD, "decode_npqtl.txt")
NPQTL <- if (file.exists(npf)) as.integer(readLines(npf)[1]) else 35559L
cat(sprintf("deCODE pQTL N = %d\n", NPQTL))

run_abf <- function(study) {
  f <- file.path(WD, paste0("coloc_input_decode_", study, ".tsv"))
  if (!file.exists(f)) { cat("MISSING", f, "\n"); return(invisible(NULL)) }
  d <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (nrow(d) < 5) { cat("too few SNPs:", study, nrow(d), "\n"); return(invisible(NULL)) }
  D1 <- list(beta = d$beta_pqtl, varbeta = d$se_pqtl^2, snp = d$snp, position = d$pos,
             type = "quant", N = NPQTL, MAF = d$maf_pqtl, sdY = 1)
  D2 <- list(beta = d$beta_gwas, varbeta = d$se_gwas^2, snp = d$snp, position = d$pos,
             type = "cc", N = as.integer(d$n_gwas[1]), MAF = d$maf_gwas, s = 0.05)
  z1 <- abs(D1$beta / sqrt(D1$varbeta)); z2 <- abs(D2$beta / sqrt(D2$varbeta))
  cat(sprintf("\n=== coloc.abf deCODE-ProteinS x %s : nsnps=%d ; pQTL max|z|=%.2f @ %d ; VTE max|z|=%.2f @ %d ===\n",
              study, nrow(d), max(z1), d$pos[which.max(z1)], max(z2), d$pos[which.max(z2)]))
  print(round(coloc.abf(D1, D2)$summary, 4))
}
run_abf("finngen"); run_abf("ukbb_neale")

cat("\n=== CONDITIONAL coloc.signals (FinnGen, 1KG-EUR LD) ===\n")
df <- file.path(WD, "cond_coloc_decode_finngen_data.tsv")
lf <- file.path(WD, "cond_coloc_decode_finngen_LD.tsv")
if (file.exists(df) && file.exists(lf)) {
  dat <- read.table(df, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  R <- as.matrix(read.table(lf, header = FALSE, sep = "\t"))
  snp <- dat$snp; pos <- dat$pos; dimnames(R) <- list(snp, snp)
  Ng <- as.integer(dat$n_gwas[1])
  D1 <- list(beta = dat$beta_pqtl, varbeta = dat$varbeta_pqtl, snp = snp, position = pos,
             type = "quant", N = NPQTL, MAF = dat$maf, sdY = 1, LD = R)
  D2 <- list(beta = dat$beta_gwas, varbeta = dat$varbeta_gwas, snp = snp, position = pos,
             type = "cc", N = Ng, MAF = dat$maf_gwas, s = 0.05, LD = R)
  z1 <- abs(D1$beta / sqrt(D1$varbeta)); z2 <- abs(D2$beta / sqrt(D2$varbeta))
  cat(sprintf("nsnps=%d ; pQTL max|z|=%.2f @ %d ; VTE max|z|=%.2f @ %d\n", length(snp),
              max(z1), pos[which.max(z1)], max(z2), pos[which.max(z2)]))
  cat("\n-- baseline coloc.abf on this SNP set --\n"); print(round(coloc.abf(D1, D2)$summary, 4))
  cat("\n-- coloc.signals method=cond --\n")
  rc <- tryCatch(coloc.signals(D1, D2, method = "cond", mode = "iterative", p12 = 5e-6),
                 error = function(e) { cat("cond ERR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(rc)) {
    s <- as.data.frame(rc$summary)
    cols <- intersect(c("hit1", "hit2", "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf"), colnames(s))
    print(s[, cols, drop = FALSE], row.names = FALSE)
    for (i in seq_len(nrow(s))) {
      h1 <- s$hit1[i]; h2 <- s$hit2[i]
      p1 <- if (h1 %in% snp) pos[match(h1, snp)] else NA
      p2 <- if (h2 %in% snp) pos[match(h2, snp)] else NA
      cat(sprintf("  row%d: pQTL-hit=%s(%s) VTE-hit=%s(%s) PP.H4=%.3f\n", i, h1, p1, h2, p2, s$PP.H4.abf[i]))
    }
  }
  cat("\n-- coloc.signals method=mask --\n")
  rmsk <- tryCatch(coloc.signals(D1, D2, method = "mask", mode = "iterative", p12 = 5e-6, r2thr = 0.01),
                   error = function(e) { cat("mask ERR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(rmsk)) {
    s <- as.data.frame(rmsk$summary)
    cols <- intersect(c("hit1", "hit2", "nsnps", "PP.H3.abf", "PP.H4.abf"), colnames(s))
    print(s[, cols, drop = FALSE], row.names = FALSE)
  }
} else cat("MISSING cond data/LD files\n")
cat("\nDONE_DECODE_COLOC\n")
