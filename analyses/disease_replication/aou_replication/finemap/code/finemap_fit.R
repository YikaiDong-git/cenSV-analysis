#!/usr/bin/env Rscript
.libPaths(c(path.expand("~/.local/lib/R"), .libPaths()))
suppressMessages(library(susieR))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[[i + 1]] }
dir <- getarg("--dir"); if (is.null(dir)) stop("need --dir")
Lmax <- as.integer(getarg("--L", "10")); coverage <- as.numeric(getarg("--coverage", "0.95"))

sh <- read.delim(file.path(dir, "shapes.tsv"), header = FALSE, col.names = c("k", "v"), stringsAsFactors = FALSE)
np <- as.integer(sh$v[sh$k == "n_present"]); V <- as.integer(sh$v[sh$k == "n_variants"])
n_eff <- as.numeric(sh$v[sh$k == "n_eff"])
meta <- read.delim(file.path(dir, "variant_meta.tsv"), stringsAsFactors = FALSE); stopifnot(nrow(meta) == V)
z_log <- as.numeric(meta$z)
R <- as.matrix(read.table(file.path(dir, "R.tsv"), header = FALSE)); stopifnot(nrow(R) == V, ncol(R) == V)
R <- (R + t(R)) / 2; diag(R) <- 1          # symmetry + unit diag ONLY (no ridge; R full-rank from residualized M)
Mres <- matrix(readBin(file.path(dir, "Mres.bin"), "double", np * V, 8, endian = "little"), np, V, byrow = TRUE)
yres <- readBin(file.path(dir, "yres.bin"), "double", np, 8, endian = "little")
L <- max(1L, min(Lmax, V - 1L))

rr <- as.numeric(cor(Mres, yres)); rr <- pmin(pmax(rr, -0.999999), 0.999999)
z_lpm <- rr * sqrt((np - 2) / (1 - rr^2))

fit_indiv <- tryCatch(susie(Mres, yres, L = L, intercept = FALSE, standardize = TRUE,
                            estimate_residual_variance = TRUE, max_iter = 1000),
                      error = function(e) { message("indiv error: ", conditionMessage(e)); NULL })
fit_rss_lpm <- tryCatch(susie_rss(z = z_lpm, R = R, n = np, L = L,
                                  estimate_residual_variance = FALSE, max_iter = 1000, check_prior = FALSE),
                        error = function(e) { message("rss_lpm error: ", conditionMessage(e)); NULL })
fit_rss_logit <- tryCatch(susie_rss(z = z_log, R = R, n = n_eff, L = L,
                                    estimate_residual_variance = FALSE, max_iter = 1000, check_prior = FALSE),
                          error = function(e) { message("rss_logit error: ", conditionMessage(e)); NULL })

set_cs <- function(f, useX) {
  if (is.null(f)) return(f)
  f$sets <- if (useX) susie_get_cs(f, X = Mres, coverage = coverage, min_abs_corr = 0.5)
            else       susie_get_cs(f, Xcorr = R, coverage = coverage, min_abs_corr = 0.5)
  f
}
fit_indiv <- set_cs(fit_indiv, TRUE); fit_rss_lpm <- set_cs(fit_rss_lpm, FALSE)
fit_rss_logit <- set_cs(fit_rss_logit, FALSE)

cs_list <- function(f) if (is.null(f) || is.null(f$sets$cs) || length(f$sets$cs) == 0) list() else f$sets$cs
pip_vec <- function(f) if (is.null(f)) rep(NA_real_, V) else as.numeric(f$pip)
cs_label <- function(csl) { lab <- rep("", V)
  if (length(csl)) for (nm in names(csl)) lab[csl[[nm]]] <- ifelse(lab[csl[[nm]]] == "", nm, paste(lab[csl[[nm]]], nm, sep = ";"))
  lab }
union_idx <- function(csl) if (length(csl)) sort(unique(unlist(csl))) else integer(0)
jac <- function(a, b) { u <- union(union_idx(a), union_idx(b))
  if (!length(u)) NA_real_ else length(intersect(union_idx(a), union_idx(b))) / length(u) }

E <- list(indiv = fit_indiv, rss_lpm = fit_rss_lpm, rss_logit = fit_rss_logit)
cs <- lapply(E, cs_list); pip <- lapply(E, pip_vec); lab <- lapply(cs, cs_label)
is_pred <- meta$is_en_predictor == 1; is_snp <- meta$variant_type == "snp"
sv_idx <- which(meta$variant_type == "sv"); snp_idx <- which(meta$variant_type == "snp")
value_gate <- function(csl) { if (!length(csl)) return(FALSE); if (length(csl) > 1) return(TRUE)
  m <- unlist(csl); any(meta$variant_type[m] == "snp" & !is_pred[m]) }
max_sv_snp_r <- if (length(sv_idx) && length(snp_idx)) max(abs(R[sv_idx, snp_idx, drop = FALSE])) else NA_real_

pips <- data.frame(variant_id = meta$variant_id, variant_type = meta$variant_type,
  is_en_predictor = as.integer(is_pred), z_logit = z_log, z_lpm = z_lpm, p = meta$p,
  pip_indiv = pip$indiv, pip_rss_lpm = pip$rss_lpm, pip_rss_logit = pip$rss_logit,
  cs_indiv = lab$indiv, cs_rss_lpm = lab$rss_lpm, cs_rss_logit = lab$rss_logit, stringsAsFactors = FALSE)
pips <- pips[order(-pmax(pips$pip_indiv, pips$pip_rss_logit, na.rm = TRUE)), ]
write.table(pips, file.path(dir, "susie_variant_pips.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cs_rows <- data.frame()
for (eng in names(E)) { csl <- cs[[eng]]; pp <- pip[[eng]]
  if (length(csl)) for (nm in names(csl)) { idx <- csl[[nm]]
    cs_rows <- rbind(cs_rows, data.frame(engine = eng, cs = nm, size = length(idx),
      n_sv = sum(meta$variant_type[idx] == "sv"), n_snp = sum(meta$variant_type[idx] == "snp"),
      n_nonpred_snp = sum(meta$variant_type[idx] == "snp" & !is_pred[idx]),
      top_variant = meta$variant_id[idx[which.max(pp[idx])]], top_type = meta$variant_type[idx[which.max(pp[idx])]],
      top_pip = max(pp[idx]), variants = paste(meta$variant_id[idx], collapse = ","), stringsAsFactors = FALSE)) } }
write.table(cs_rows, file.path(dir, "susie_credible_sets.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

top_sv_pip <- function(pp) if (length(sv_idx)) max(pp[sv_idx]) else NA_real_
sum_sv_pip <- function(pp) if (length(sv_idx)) sum(pp[sv_idx]) else NA_real_
summ <- data.frame(n_variants = V, n_sv = length(sv_idx), n_snp = length(snp_idx), n_present = np, n_eff = n_eff, L = L,
  n_cs_indiv = length(cs$indiv), n_cs_rss_lpm = length(cs$rss_lpm), n_cs_rss_logit = length(cs$rss_logit),
  jaccard_indiv_vs_rsslpm = jac(cs$indiv, cs$rss_lpm),
  pipcor_indiv_vs_rsslpm = if (all(is.finite(pip$indiv)) && all(is.finite(pip$rss_lpm))) cor(pip$indiv, pip$rss_lpm) else NA_real_,
  jaccard_indiv_vs_rsslogit = jac(cs$indiv, cs$rss_logit),
  pipcor_indiv_vs_rsslogit = if (all(is.finite(pip$indiv)) && all(is.finite(pip$rss_logit))) cor(pip$indiv, pip$rss_logit) else NA_real_,
  value_indiv = value_gate(cs$indiv), value_rss_lpm = value_gate(cs$rss_lpm), value_rss_logit = value_gate(cs$rss_logit),
  sum_sv_pip_indiv = sum_sv_pip(pip$indiv), sum_sv_pip_rss_lpm = sum_sv_pip(pip$rss_lpm), sum_sv_pip_rss_logit = sum_sv_pip(pip$rss_logit),
  top_sv_pip_indiv = top_sv_pip(pip$indiv), top_sv_pip_rss_logit = top_sv_pip(pip$rss_logit),
  max_abs_sv_snp_r = max_sv_snp_r,
  converged_indiv = !is.null(fit_indiv) && isTRUE(fit_indiv$converged),
  converged_rss_lpm = !is.null(fit_rss_lpm) && isTRUE(fit_rss_lpm$converged),
  converged_rss_logit = !is.null(fit_rss_logit) && isTRUE(fit_rss_logit$converged), stringsAsFactors = FALSE)
write.table(summ, file.path(dir, "susie_fit_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("FIT %s | V=%d (sv=%d snp=%d) | nCS indiv=%d rss_lpm=%d rss_logit=%d | QC indiv~rss_lpm: Jac=%s PIPcor=%.4f | "
    , basename(dir), V, length(sv_idx), length(snp_idx), length(cs$indiv), length(cs$rss_lpm), length(cs$rss_logit),
    ifelse(is.na(summ$jaccard_indiv_vs_rsslpm), "NA", sprintf("%.3f", summ$jaccard_indiv_vs_rsslpm)), summ$pipcor_indiv_vs_rsslpm))
cat(sprintf("value indiv=%s | topSVpip indiv=%.3f logit=%.3f | max|r(SV,SNP)|=%.4f\n",
    summ$value_indiv, top_sv_pip(pip$indiv), top_sv_pip(pip$rss_logit), ifelse(is.na(max_sv_snp_r), NA, max_sv_snp_r)))
