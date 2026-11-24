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
tsupport <- as.integer(b[["tsupport_local"]]) + 1L
sv_row <- p + tcol + 1L; z_sv <- z_aug[sv_row]
cat(sprintf("p=%d uncond z_SV=%.3f n_med=%.0f n_support=%d\n", p, z_sv, n_med, length(tsupport)))

gene_of <- function(pp) {
  if (pp >= 93873038 && pp <= 93980003) return("PROS1")
  if (pp >= 93980128 && pp <= 94055678) return("ARL13B")
  if (pp >= 94062980 && pp <= 94131832) return("NSUN3")
  if (pp >= 89107522 && pp <= 89482437) return("EPHA3")
  if (pp <  88500000) return("prox(87)")
  if (pp >= 93800000) return("PROS1nbhd")
  if (pp >= 93000000) return("distPROS1nbhd")
  return("interior")
}
cond_single_aug <- function(j) { rj <- R_aug[sv_row, j]; (z_sv - rj * z_aug[j]) / sqrt(max(1 - rj^2, 1e-6)) }
purity_of <- function(idx) { if (length(idx) <= 1) return(1.0)
  sub <- abs(Roo[idx, idx]); min(sub[lower.tri(sub)]) }

fit <- susie_rss(z = zoo, R = Roo, L = 10L, max_iter = 200L, coverage = 0.95,
                 estimate_residual_variance = FALSE, check_prior = FALSE)
A <- fit$alpha                      # L x p posterior assignment per single effect
L <- nrow(A)

cat("\n===== RAW per-effect structure (BEFORE purity filter) =====\n")
cat("effect: MAP variant (top alpha) | its 95% inclusion set size+purity | link to sv40\n")
for (l in 1:L) {
  a <- A[l, ]; top <- which.max(a)
  ord <- order(-a); cum <- cumsum(a[ord]); k <- which(cum >= 0.95)[1]; csl <- ord[1:k]
  cat(sprintf("  eff%-2d top pos=%d gene=%-13s alpha=%.3f z=%+.2f | CS95 size=%-4d purity=%.3f | r(sv40,top)=%+.3f  cond z_SV %.2f->%.2f\n",
              l, snp_pos[top], gene_of(snp_pos[top]), a[top], zoo[top], length(csl),
              purity_of(csl), R_aug[sv_row, top], z_sv, cond_single_aug(top)))
}

cat("\n===== credible sets at RELAXED purity (min_abs_corr) =====\n")
for (mac in c(0.5, 0.1, 0.0)) {
  cs <- susieR::susie_get_cs(fit, Xcorr = Roo, coverage = 0.95, min_abs_corr = mac)
  pip <- susieR::susie_get_pip(fit)
  cat(sprintf("\n-- min_abs_corr=%.2f -> %d CS --\n", mac, length(cs$cs)))
  if (length(cs$cs) > 0) for (k in seq_along(cs$cs)) {
    idx <- cs$cs[[k]]; ld <- idx[which.max(pip[idx])]
    cat(sprintf("   CS%-2d %-7s size=%-4d lead pos=%d gene=%-13s z=%+.2f purity=%.3f | r(sv40,lead)=%+.3f cond z_SV %.2f->%.2f overlap_sv40supp=%d\n",
                k, names(cs$cs)[k], length(idx), snp_pos[ld], gene_of(snp_pos[ld]), zoo[ld],
                purity_of(idx), R_aug[sv_row, ld], z_sv, cond_single_aug(ld), sum(idx %in% tsupport)))
  }
}

cat("\n===== PROS1-neighborhood (93.8-94.4 Mb) direct look =====\n")
nbhd <- which(snp_pos >= 93800000 & snp_pos < 94400000)
if (length(nbhd) > 0) {
  rsv <- abs(R_aug[sv_row, nbhd])
  cat(sprintf("  #SNPs in PROS1-nbhd support=%d ; with |r(sv40)|>0.3 : %d ; >0.2 : %d ; max|r|=%.3f\n",
              length(nbhd), sum(rsv > 0.3), sum(rsv > 0.2), max(rsv)))
  ord <- nbhd[order(-rsv[order(match(nbhd, nbhd))])]   # keep simple: sort nbhd by |r|
  ord <- nbhd[order(-abs(R_aug[sv_row, nbhd]))]
  cat("  top 6 PROS1-nbhd SNPs by |r(sv40)|:\n")
  for (j in head(ord, 6)) cat(sprintf("    pos=%d gene=%-13s z=%+.2f r(sv40)=%+.3f cond z_SV->%.2f\n",
              snp_pos[j], gene_of(snp_pos[j]), zoo[j], R_aug[sv_row, j], cond_single_aug(j)))
  cond_joint_aug <- function(C){ Rcc<-R_aug[C,C,drop=FALSE]; Rsc<-R_aug[sv_row,C,drop=FALSE]
    s<-svd(Rcc); di<-ifelse(s$d>1e-10*max(s$d),1/s$d,0); inv<-s$v%*%(di*t(s$u))
    num<-z_sv-as.numeric(Rsc%*%inv%*%z_aug[C]); v<-max(1-as.numeric(Rsc%*%inv%*%t(Rsc)),1e-6); num/sqrt(v) }
  cat("  sequential JOINT conditioning on the top-k |r| PROS1-nbhd SNPs:\n")
  for (k in c(1,2,3,5,10)) { kk<-min(k,length(ord)); cat(sprintf("    k=%-2d : z_SV %.2f -> %.2f\n", kk, z_sv, cond_joint_aug(ord[1:kk]))) }
}
cat("\nDONE_CS_PURITY_PROBE\n")
