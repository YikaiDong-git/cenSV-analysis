#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(reticulate) })
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

VENDOR <- "${CENSV_ROOT}/analyses/sv_catalog/code/q8_validation/qc/R/vendor_susieR/susie_rss_lambda_0_14_2.R"
source(VENDOR)
if (as.character(packageVersion("susieR")) != "0.14.2") stop("need susieR 0.14.2")
np <- import("numpy")

args <- commandArgs(trailingOnly = TRUE)
manifest <- args[1]; out_tsv <- args[2]
smoke <- if (length(args) >= 3 && nzchar(args[3])) strsplit(args[3], ",")[[1]] else NULL

M <- read.table(manifest, sep = "\t", stringsAsFactors = FALSE, quote = "")
colnames(M) <- c("jw", "mode", "arm", "trait", "region", "hint", "col")

con <- file(out_tsv, "w")
writeLines(paste("region", "trait", "arm", "mode", "sv_locus", "sv_index", "cv_r2", "sv_pip",
                 "svSNP_r", "z_SV", "in_cs", "p_snps", "G", "converged", "niter", "wall_s",
                 sep = "\t"), con)

n_ok <- 0
for (i in seq_len(nrow(M))) {
  jw <- M$jw[i]
  if (!file.exists(jw)) next
  b <- np$load(jw, allow_pickle = FALSE)
  R_aug <- b[["R_aug"]]; z_aug <- as.numeric(b[["z_aug"]])
  sv_locus <- as.character(b[["sv_locus"]]); sv_index <- as.integer(b[["sv_index"]])
  cv <- as.numeric(b[["sv_cv_r2"]])
  G <- length(sv_locus); p <- length(z_aug) - G
  if (p < 3 || G < 1) next
  t0 <- Sys.time()
  fit <- tryCatch(
    susie_rss(z = z_aug, R = R_aug, L = 10L, max_iter = 200L,
              estimate_residual_variance = FALSE, check_prior = FALSE),
    error = function(e) NULL)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(fit)) next
  n_ok <- n_ok + 1
  pip <- fit$pip; conv <- isTRUE(fit$converged); niter <- fit$niter %||% NA
  cs <- tryCatch(susie_get_cs(fit, coverage = 0.95, min_abs_corr = 0.5)$cs, error = function(e) list())
  in_cs_vars <- unlist(cs)
  for (g in seq_len(G)) {
    loc <- sv_locus[g]
    if (M$mode[i] == "allwindow" && !is.null(smoke) && !(loc %in% smoke)) next
    col <- p + g
    svsnp_r <- max(abs(R_aug[col, 1:p]))
    writeLines(paste(M$region[i], M$trait[i], M$arm[i], M$mode[i], loc, sv_index[g],
                     sprintf("%.4f", cv[g]), sprintf("%.5f", pip[col]), sprintf("%.4f", svsnp_r),
                     sprintf("%.4f", z_aug[col]), as.integer(col %in% in_cs_vars),
                     p, G, conv, niter, sprintf("%.2f", wall), sep = "\t"), con)
  }
}
close(con)
cat(sprintf("[run_susie_batch] %d/%d bundles fit -> %s\n", n_ok, nrow(M), out_tsv))
