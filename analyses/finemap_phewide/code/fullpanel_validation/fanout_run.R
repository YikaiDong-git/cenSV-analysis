#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(reticulate); library(susieR) })
args <- commandArgs(TRUE)
REGION <- args[1]; MANIFEST <- args[2]; OUTDIR <- args[3]; CODE <- args[4]
CSTART <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else NA
CEND   <- if (length(args) >= 6 && nzchar(args[6])) as.integer(args[6]) else NA
CHUNK  <- if (length(args) >= 7 && nzchar(args[7])) args[7] else "all"
k <- import_from_path("fanout_kernel", CODE)

man <- read.table(MANIFEST, header = TRUE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
rows <- man[man$region == REGION, ]
if (nrow(rows) == 0) { cat(sprintf("[fanout] no cells for %s\n", REGION)); quit(status = 0) }
if (!is.na(CSTART)) { rows <- rows[CSTART:min(CEND, nrow(rows)), , drop = FALSE]; rows <- rows[!is.na(rows$region), , drop = FALSE] }
chrom <- rows$chrom[1]; rs <- rows$region_start[1]; re <- rows$region_end[1]
cg <- sub("chr", "", chrom)

t0 <- Sys.time()
p_full <- k$load_region(chrom, rs, re)
cat(sprintf("[fanout] %s loaded: p_full=%d cells=%d (%.0fs)\n", REGION, p_full, nrow(rows),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

outf <- file.path(OUTDIR, paste0(REGION, "_c", CHUNK, "_fanout_diag.tsv"))
con <- file(outf, "w")
cat("region\tsv_locus\tsv_index\ttrait\tarm\tcv\tsvSNP_r_union\tpip_matched_union\tpip_impute_union\tp_snp\tsv_pip\tmax_snp_pip\tsv_in_cs\tn_cs\tconverged\tz_SV\tsecs\n", file = con)

for (i in seq_len(nrow(rows))) {
  r <- rows[i, ]
  for (arm in c("matched", "impute")) {
    base <- sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", REGION, r$sv_locus, r$sv_index, r$trait, arm,
                    r$cv, r$svSNP_r, r$pip_matched, r$pip_impute)
    cell <- tryCatch(k$build_cell(r$sv_index, r$sv_locus, r$trait, r$bgz, cg, rs, re, arm, 12000L),
                     error = function(e) { cat(sprintf("[err build] %s %s %s: %s\n", r$sv_locus, r$trait, arm, conditionMessage(e))); NULL })
    if (is.null(cell)) { cat(sprintf("%s\tSKIP\tNA\tNA\tNA\tNA\tNA\tNA\n", base), file = con); flush(con); next }
    n <- cell$n
    R <- matrix(as.numeric(cell$R_aug), n, n)   # symmetric: byrow irrelevant
    z <- as.numeric(cell$z_aug)
    tt <- Sys.time()
    fit <- tryCatch({
      setTimeLimit(elapsed = 180, transient = TRUE)
      on.exit(setTimeLimit(cpu = Inf, elapsed = Inf), add = TRUE)
      susie_rss(z = z, R = R, L = 10, max_iter = 100,
                estimate_residual_variance = FALSE, check_prior = FALSE)
    }, error = function(e) { cat(sprintf("[err susie] %s %s %s: %s\n", r$sv_locus, r$trait, arm, conditionMessage(e))); NULL })
    if (is.null(fit)) { cat(sprintf("%s\t%d\tSUSIE_ERR\tNA\tNA\tNA\tNA\t%.2f\tNA\n", base, cell$p_snp, cell$z_SV), file = con); flush(con); next }
    dt <- as.numeric(difftime(Sys.time(), tt, units = "secs"))
    pip <- fit$pip; sv <- pip[n]
    cs <- tryCatch(susie_get_cs(fit, coverage = 0.95, min_abs_corr = 0.5)$cs, error = function(e) list())
    incs <- any(vapply(cs, function(s) n %in% s, logical(1)))
    cat(sprintf("%s\t%d\t%.5f\t%.4f\t%s\t%d\t%s\t%.2f\t%.0f\n", base, cell$p_snp, sv,
                max(pip[1:(n - 1)]), incs, length(cs), fit$converged, cell$z_SV, dt), file = con)
    flush(con)
  }
}
close(con)
cat(sprintf("[fanout] %s DONE (%.0fmin)\n", REGION, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
