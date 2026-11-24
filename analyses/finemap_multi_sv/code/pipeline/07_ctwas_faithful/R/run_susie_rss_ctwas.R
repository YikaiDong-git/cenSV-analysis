suppressPackageStartupMessages({
  library(reticulate)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || (length(a) == 1L && is.na(a))) b else a

if (!requireNamespace("susieR", quietly = TRUE)) {
  stop("susieR package not installed")
}
.SUSIE_PINNED_VERSION <- "0.14.2"
sv_pkg_ver <- as.character(packageVersion("susieR"))
if (sv_pkg_ver != .SUSIE_PINNED_VERSION) {
  stop(sprintf("susieR version mismatch: have %s want %s", sv_pkg_ver, .SUSIE_PINNED_VERSION))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: Rscript run_susie_rss_ctwas.R <jw.npz> <out_dir> <pheno> [<region_id>]")
jw_path   <- args[1]
out_dir   <- args[2]
pheno     <- args[3]
region_id <- if (length(args) >= 4) args[4] else sub("__.*$", "", basename(jw_path))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

np <- import("numpy")
bundle <- np$load(jw_path, allow_pickle = FALSE)

R_aug    <- bundle[["R_aug"]]
z_aug    <- as.numeric(bundle[["z_aug"]])
sv_locus    <- as.character(bundle[["sv_locus"]])
sv_index    <- as.integer(bundle[["sv_index"]])
sv_category <- as.character(bundle[["sv_category"]])
sv_strategy <- as.character(bundle[["sv_strategy"]])
sv_arm      <- as.character(bundle[["sv_arm"]])
sv_cv_r2    <- as.numeric(bundle[["sv_cv_r2"]])
G <- length(sv_locus)
p <- length(z_aug) - G
N_GWAS <- length(z_aug)  # cTWAS convention: n = length(z)

cat(sprintf("[ctwas-susie] region=%s pheno=%s p=%d G=%d susieR=%s L=10 max_iter=100 n=%d\n",
            region_id, pheno, p, G, sv_pkg_ver, N_GWAS))

stopifnot(nrow(R_aug) == ncol(R_aug),
          nrow(R_aug) == p + G,
          length(z_aug) == p + G,
          length(sv_locus) == G)

t0 <- Sys.time()
fit <- susieR::susie_rss(
  z = z_aug,
  R = R_aug,
  n = N_GWAS,                                # cTWAS-faithful: n=length(z)
  L = 10L,                                    # cTWAS default (vs Round-2 fix L=5)
  max_iter = 100L,                            # cTWAS vendored default (vs Round-2 fix max_iter=500)
  coverage = 0.95,
  estimate_residual_variance = FALSE,
  estimate_prior_variance = FALSE,
  prior_variance = 50,
  check_prior = FALSE
)
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

pip <- susieR::susie_get_pip(fit)
stopifnot(length(pip) == p + G)
pip_snp <- pip[seq_len(p)]
pip_sv  <- pip[(p + 1L):(p + G)]

n_iter    <- fit$niter %||% NA_integer_
converged <- fit$converged %||% NA
sigma2    <- fit$sigma2 %||% NA_real_

diag_tsv <- data.frame(
  region_id   = region_id,
  pheno       = pheno,
  sv_locus    = sv_locus,
  sv_index    = sv_index,
  sv_category = sv_category,
  sv_strategy = sv_strategy,
  sv_arm      = sv_arm,
  cv_r2       = sv_cv_r2,
  sv_pip      = pip_sv,
  p_snps      = p,
  G_svs       = G,
  n_iter      = as.integer(n_iter),
  converged   = as.character(converged),
  sigma2      = sigma2,
  susieR_version = sv_pkg_ver,
  L = 10L, max_iter = 100L, check_prior = FALSE,
  n_passed = N_GWAS,
  prior_variance = 50,
  estimate_prior_variance = FALSE,
  variant = "ctwas_faithful",
  r_aug_mode = "canonical_diag1_normalized",
  wall_seconds = dt,
  stringsAsFactors = FALSE
)

out_diag <- file.path(out_dir, sprintf("%s__%s__diag.tsv", region_id, pheno))
tmp <- paste0(out_diag, sprintf(".tmp.%d", Sys.getpid()))
write.table(diag_tsv, tmp, sep = "\t", quote = FALSE, row.names = FALSE)
file.rename(tmp, out_diag)

out_snp <- file.path(out_dir, sprintf("%s__%s__pip_snp.tsv", region_id, pheno))
snp_id    <- as.character(bundle[["snp_id"]])
snp_chrom <- as.character(bundle[["snp_chrom"]])
snp_pos   <- as.integer(bundle[["snp_pos"]])
snp_df <- data.frame(snp_id = snp_id, snp_chrom = snp_chrom, snp_pos = snp_pos,
                     snp_pip = pip_snp, stringsAsFactors = FALSE)
tmp2 <- paste0(out_snp, sprintf(".tmp.%d", Sys.getpid()))
write.table(snp_df, tmp2, sep = "\t", quote = FALSE, row.names = FALSE)
file.rename(tmp2, out_snp)

cat(sprintf("[ctwas-susie] DONE p=%d G=%d wall=%.1fs sv_pip max=%.4f n_sv>=0.1=%d conv=%s niter=%d\n",
            p, G, dt, max(pip_sv), sum(pip_sv >= 0.1), as.character(converged), as.integer(n_iter)))
