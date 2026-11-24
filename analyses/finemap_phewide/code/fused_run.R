#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(reticulate) })
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
VENDOR <- "${CENSV_ROOT}/analyses/sv_catalog/code/q8_validation/qc/R/vendor_susieR/susie_rss_lambda_0_14_2.R"
CODE   <- "${CENSV_ROOT}/analyses/finemap_phewide/code"
source(VENDOR)
if (as.character(packageVersion("susieR")) != "0.14.2") stop("need susieR 0.14.2")
k <- import_from_path("fused_kernel", CODE)

a <- commandArgs(trailingOnly = TRUE)
bundle <- a[1]; zmat <- a[2]; mode <- a[3]; arm <- a[4]
cstart <- as.integer(a[5]); cend <- as.integer(a[6]); region <- a[7]; out <- a[8]

dims <- k$load(bundle, zmat); Tn <- dims[[3]]
traits <- k$traits(); nobs <- unlist(k$nobs())
con <- file(out, "w")
writeLines(paste("region", "trait", "arm", "mode", "sv_locus", "sv_index", "cv_r2", "sv_pip",
                 "svSNP_r", "z_SV", "in_cs", "p_snps", "m_sv", "s_rss", "converged", "niter", "wall_s",
                 sep = "\t"), con)

process_fit <- function(R_aug, z_aug, sv_meta, trait) {
  m <- length(sv_meta); if (m < 1 || is.null(R_aug)) return(invisible())   # defensive: never fit an empty panel
  pj <- nrow(R_aug) - m
  t0 <- Sys.time()
  fit <- tryCatch(susie_rss(z = z_aug, R = R_aug, L = 10L, max_iter = 200L,
                            estimate_residual_variance = FALSE, check_prior = FALSE),
                  error = function(e) NULL)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(fit)) return(invisible())
  cs <- tryCatch(susie_get_cs(fit, coverage = 0.95, min_abs_corr = 0.5)$cs, error = function(e) list())
  incs <- unlist(cs)
  svcols <- (pj + 1):(pj + m)
  s_fit <- if (isTRUE(any(abs(z_aug[svcols]) > 4, na.rm = TRUE))) tryCatch(estimate_s_rss(z = z_aug, R = R_aug), error = function(e) NA_real_) else NA_real_
  conv <- isTRUE(fit$converged); niter <- fit$niter %||% NA
  for (g in seq_len(m)) {
    col <- pj + g; meta <- sv_meta[[g]]
    writeLines(paste(region, trait, arm, mode, meta[[1]], meta[[2]], sprintf("%.4f", meta[[3]]),
                     sprintf("%.5f", fit$pip[col]), sprintf("%.4f", max(abs(R_aug[col, 1:pj]))),
                     sprintf("%.4f", z_aug[col]), as.integer(col %in% incs), pj, m,
                     ifelse(is.na(s_fit), "NA", sprintf("%.4f", s_fit)),
                     conv, niter, sprintf("%.2f", wall), sep = "\t"), con)
  }
}

if (mode == "hitbundle") {
  for (j in cstart:(cend - 1)) {
    meta <- list(k$sv_info(as.integer(j)))
    if (arm == "impute") {
      R_aug <- k$hb_impute_raug(as.integer(j)); if (is.null(R_aug)) next
      for (t in 0:(Tn - 1)) {
        if (nobs[t + 1] < 5) next
        process_fit(R_aug, as.numeric(k$hb_impute_zaug(as.integer(t))), meta, traits[t + 1])
      }
    } else {
      for (t in 0:(Tn - 1)) {
        if (nobs[t + 1] < 5) next
        rr <- k$hb_matched(as.integer(j), as.integer(t)); if (is.null(rr)) next
        process_fit(rr[[1]], as.numeric(rr[[2]]), meta, traits[t + 1])
      }
    }
  }
} else {                                              # all_window
  if (arm == "impute") {
    R_aug <- k$aw_impute_raug()
    if (!is.null(R_aug)) {
      meta <- lapply(unlist(k$aw_impute_keep()), function(jj) k$sv_info(as.integer(jj)))
      for (t in 0:(Tn - 1)) {
        if (nobs[t + 1] < 5) next
        process_fit(R_aug, as.numeric(k$aw_impute_zaug(as.integer(t))), meta, traits[t + 1])
      }
    }
  } else {
    for (t in 0:(Tn - 1)) {
      if (nobs[t + 1] < 5) next
      rr <- k$aw_matched(as.integer(t)); if (is.null(rr)) next
      meta <- lapply(unlist(rr[[3]]), function(jj) k$sv_info(as.integer(jj)))
      process_fit(rr[[1]], as.numeric(rr[[2]]), meta, traits[t + 1])
    }
  }
}
close(con)
cat(sprintf("[fused_run] %s %s %s cols[%d,%d) -> %s\n", region, mode, arm, cstart, cend, out))
