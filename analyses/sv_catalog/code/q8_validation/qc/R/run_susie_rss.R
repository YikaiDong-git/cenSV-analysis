suppressPackageStartupMessages({
  library(data.table)
  library(reticulate)
})

.cargs <- commandArgs(trailingOnly = FALSE)
.fileflag <- .cargs[grep("^--file=", .cargs)]
.HERE <- Sys.getenv("V2_CANONICAL_QC_R_DIR", unset = "")
if (!nzchar(.HERE)) {
  if (length(.fileflag) == 1L) {
    .HERE <- dirname(normalizePath(sub("^--file=", "", .fileflag), mustWork = FALSE))
  } else {
    .HERE <- file.path(getwd(), "R")
  }
}
Sys.setenv(V2_CANONICAL_QC_R_DIR = .HERE)
source(file.path(.HERE, "fit_susie_rss_centromeric.R"))

args <- commandArgs(trailingOnly = TRUE)
parse_flag <- function(name) {
  i <- which(args == paste0("--", name))
  if (length(i) != 1) stop("flag --", name, " missing or duplicated")
  args[i + 1L]
}
BUNDLE   <- parse_flag("bundle-dir")
PHENO    <- parse_flag("pheno")
SV_LOCUS <- parse_flag("sv-locus")
ARM      <- parse_flag("arm")
CHROM    <- parse_flag("chrom")
OUT_DIR  <- parse_flag("out-dir")

if (!startsWith(CHROM, "chr")) CHROM <- paste0("chr", CHROM)

np <- reticulate::import("numpy")
z_aug <- as.numeric(np$load(file.path(BUNDLE, "z_aug.npy")))
R_aug <- np$load(file.path(BUNDLE, "R_aug.npy"))
R_aug <- as.matrix(R_aug)

meta_dt <- fread(file.path(BUNDLE, "meta.tsv"), sep = "\t", header = TRUE)
meta <- setNames(as.list(as.character(meta_dt$value)), meta_dt$key)
N_GWAS <- as.integer(meta$N_GWAS)
stopifnot(N_GWAS > 0)

cat(sprintf("[load] bundle=%s p=%d N_GWAS=%d Z_SV=%s wRw=%s eig_pre=%s\n",
            basename(BUNDLE), length(z_aug), N_GWAS, meta$Z_SV, meta$wRw,
            meta$eig_pre))

t0 <- Sys.time()
res <- fit_susie_rss_centromeric(
  z = z_aug, R = R_aug,
  n_gwas = N_GWAS,           # for diagnostics only
  bundle_meta = meta
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

fit <- res$fit
diag_extra <- res$diag

pip_vec <- fit$pip
sv_pip <- as.numeric(pip_vec[length(pip_vec)])
top_snp_pip <- if (length(pip_vec) > 1L) max(pip_vec[1:(length(pip_vec) - 1L)]) else NA_real_
sum_pip <- sum(pip_vec)

cat(sprintf("[fit] elapsed=%.2fs niter=%d converged=%s sv_pip=%.6f top_snp_pip=%.6f sum_pip=%.3f\n",
            elapsed, diag_extra$niter, diag_extra$converged,
            sv_pip, top_snp_pip, sum_pip))
cat(sprintf("[qc] s_rss=%.4f tier=%s kriging_outliers=%d kriging_logBF_SV=%.3f z_colspace_residual=%.3e\n",
            diag_extra$s_rss %||% NA, diag_extra$s_rss_tier %||% NA,
            diag_extra$kriging_outlier_count %||% NA_integer_,
            diag_extra$kriging_logBF_SV %||% NA,
            diag_extra$z_colspace_residual %||% NA))

diag_table <- data.table(
  phenotype = PHENO,
  sv_locus = SV_LOCUS,
  arm = ARM,
  chrom = CHROM,
  L = 10L,
  bundle_dir = basename(BUNDLE),
  n_LD_variants = length(z_aug) - 1L,
  Z_SV = as.numeric(meta$Z_SV),
  N_GWAS = N_GWAS,
  wRw = as.numeric(meta$wRw),
  eig_pre = as.numeric(meta$eig_pre),
  eig_post = as.numeric(meta$eig_post),
  stabilized = meta$stabilized,
  r_diagonal_mode = meta$r_diagonal_mode,
  r_diagonal_multiplier = as.numeric(meta$r_diagonal_multiplier),
  sv_pip = sv_pip,
  top_snp_pip = top_snp_pip,
  sum_pip = sum_pip,
  niter = diag_extra$niter,
  converged = diag_extra$converged,
  near_cap_flag = diag_extra$near_cap_flag,
  warn_unconverged = diag_extra$warn_unconverged,
  s_rss = diag_extra$s_rss,
  s_rss_tier = diag_extra$s_rss_tier,
  kriging_outlier_count = diag_extra$kriging_outlier_count,
  kriging_logBF_SV = diag_extra$kriging_logBF_SV,
  z_colspace_residual = diag_extra$z_colspace_residual,
  z_colspace_rank = diag_extra$z_colspace_rank,
  check_input_max_abs_z = diag_extra$check_input_max_abs_z,
  check_input_max_abs_R_offdiag = diag_extra$check_input_max_abs_R_offdiag,
  check_input_max_diag_dev = diag_extra$check_input_max_diag_dev,
  check_input_symmetric = diag_extra$check_input_symmetric,
  susieR_version = diag_extra$susieR_version,
  sessionInfo_hash = diag_extra$sessionInfo_hash,
  R_version = diag_extra$R_version,
  wrapper_id = diag_extra$wrapper_id,
  elapsed_sec = elapsed
)

out_dir <- file.path(OUT_DIR, PHENO, CHROM)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_diag <- file.path(out_dir, paste0(SV_LOCUS, "__", ARM, "__diag.tsv"))
tmp_diag <- paste0(out_diag, ".tmp.", Sys.getpid())
fwrite(diag_table, tmp_diag, sep = "\t")
file.rename(tmp_diag, out_diag)
cat(sprintf("[OUT] wrote %s\n", out_diag))
