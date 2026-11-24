suppressPackageStartupMessages({ library(reticulate); library(susieR) })
P    <- "${CENSV_WORKSPACE}"
BASE <- file.path(P, "Project1_Centromere/analyses/finemap_multi_sv")
VEND <- file.path(BASE, "../sv_catalog/code/q8_validation/qc/R/vendor_susieR/susie_rss_lambda_0_14_2.R")
OUTDIR <- file.path(BASE, "validation/ssimp_impg/work/diag_chr3")
source(VEND)
cat(sprintf("susieR pkg version: %s\n", as.character(packageVersion("susieR"))))

np <- import("numpy")
b  <- np$load(file.path(OUTDIR, "diag_chr3_matched.npz"), allow_pickle = TRUE)
R_aug <- b[["R_aug"]]; z_aug <- as.numeric(b[["z_aug"]])
Roo <- b[["Roo"]]; zoo <- as.numeric(b[["zoo"]]); Noo <- as.numeric(b[["Noo"]])
n_med <- as.numeric(median(Noo, na.rm = TRUE))
p <- as.integer(b[["p"]]); G <- as.integer(b[["G"]]); tcol <- as.integer(b[["tcol"]])
snp_pos <- as.integer(b[["snp_pos_obs"]])
tsupport <- as.integer(b[["tsupport_local"]]) + 1L      # py 0-based -> R 1-based (into SNP block)
sv_row <- p + tcol + 1L                                 # SV's row in augmented (1-based)
z_sv <- z_aug[sv_row]
cat(sprintf("p=%d G=%d tcol=%d sv_row=%d  uncond z_SV=%.3f  n_support=%d\n",
            p, G, tcol, sv_row, z_sv, length(tsupport)))

gene_of <- function(pp) {
  if (pp >= 93873038 && pp <= 93980003) return("PROS1")
  if (pp >= 93980128 && pp <= 94055678) return("ARL13B")
  if (pp >= 94062980 && pp <= 94131832) return("NSUN3")
  if (pp >= 89107522 && pp <= 89482437) return("EPHA3")
  if (pp <  88500000) return("prox(87)")
  if (pp >= 93000000) return("distPROS1nbhd")
  return("interior")
}
ginv2 <- function(M, tol = 1e-10) {
  s <- svd(M); d <- s$d; di <- ifelse(d > tol * max(d), 1 / d, 0)
  s$v %*% (di * t(s$u))
}
cond_single_aug <- function(j) {
  rj <- R_aug[sv_row, j]; (z_sv - rj * z_aug[j]) / sqrt(max(1 - rj^2, 1e-6))
}
cond_joint_aug <- function(C) {
  Rcc <- R_aug[C, C, drop = FALSE]; Rsc <- R_aug[sv_row, C, drop = FALSE]
  inv <- ginv2(Rcc)
  num <- z_sv - as.numeric(Rsc %*% inv %*% z_aug[C])
  v <- 1 - as.numeric(Rsc %*% inv %*% t(Rsc)); v <- max(v, 1e-6)
  num / sqrt(v)
}

analyze <- function(n_arg, tag) {
  cat(sprintf("\n############## %s ##############\n", tag))
  fit <- if (is.null(n_arg))
    susie_rss(z = zoo, R = Roo, L = 10L, max_iter = 200L, coverage = 0.95,
              estimate_residual_variance = FALSE, check_prior = FALSE)
  else
    susie_rss(z = zoo, R = Roo, n = n_arg, L = 10L, max_iter = 200L, coverage = 0.95,
              estimate_residual_variance = FALSE, check_prior = FALSE)
  pip <- susieR::susie_get_pip(fit)
  cs  <- susieR::susie_get_cs(fit, Xcorr = Roo, coverage = 0.95, min_abs_corr = 0.5)
  csets <- cs$cs
  cat(sprintf("SNP-only SuSiE-RSS: %d credible set(s) (L=10, coverage 0.95, purity>=0.5)\n",
              length(csets)))
  leads <- integer(0)
  if (length(csets) > 0) for (k in seq_along(csets)) {
    idx <- csets[[k]]; ld <- idx[which.max(pip[idx])]; leads <- c(leads, ld)
    cat(sprintf("  CS%-2d %-7s size=%-4d lead pos=%d z=%+.2f gene=%-13s PIP_lead=%.3f  overlap_sv40support=%d\n",
                k, names(csets)[k], length(idx), snp_pos[ld], zoo[ld], gene_of(snp_pos[ld]),
                pip[ld], sum(idx %in% tsupport)))
  }
  if (length(leads) > 0) {
    cat("-- condition z_SV on each CS lead (single) --\n")
    for (k in seq_along(leads)) {
      j <- leads[k]
      cat(sprintf("  cond CS%-2d lead pos=%d gene=%-13s : z_SV %.3f -> %6.3f | r(sv40,lead)=%+.3f\n",
                  k, snp_pos[j], gene_of(snp_pos[j]), z_sv, cond_single_aug(j), R_aug[sv_row, j]))
    }
    cat(sprintf("  JOINT on all %d CS leads        : z_SV %.3f -> %6.3f\n",
                length(leads), z_sv, cond_joint_aug(leads)))
    uni <- sort(unique(unlist(csets)))
    cat(sprintf("  JOINT on FULL CS union (%d SNPs) : z_SV %.3f -> %6.3f\n",
                length(uni), z_sv, cond_joint_aug(uni)))
    ext <- leads[!(leads %in% tsupport)]
    cat(sprintf("  JOINT on CS leads NOT in sv40 support (%d) : z_SV %.3f -> %6.3f  (signal from OTHER loci)\n",
                length(ext), z_sv, if (length(ext) > 0) cond_joint_aug(ext) else z_sv))
  }
}
cat(sprintf("\nn_med (median GWAS N over observed support) = %.0f\n", n_med))
analyze(NULL,  "FIT A: production config (no n; assumes n=Inf, conservative)")
analyze(n_med, sprintf("FIT B: with n = %.0f (more power to split signals)", n_med))

cat(sprintf("\n== [sanity, tautological] JOINT on sv40's own %d support SNPs : z_SV %.3f -> %6.3f (expect ~0)\n",
            length(tsupport), z_sv, cond_joint_aug(tsupport)))
cat("\nDONE_CS_CONDITION\n")
