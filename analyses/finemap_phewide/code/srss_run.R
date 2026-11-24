#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(reticulate); library(susieR) })
if (as.character(packageVersion("susieR")) != "0.14.2") stop("need susieR 0.14.2")
np <- import("numpy")
args <- commandArgs(trailingOnly = TRUE); manifest <- args[1]; out <- args[2]
M <- read.table(manifest, sep = "\t", stringsAsFactors = FALSE, quote = "")
colnames(M) <- c("npz", "phys", "trait", "arm", "region", "reason", "p")
con <- file(out, "w")
writeLines(paste("phys", "trait", "arm", "region", "reason", "p", "s_rss", "wall_s", sep = "\t"), con)
for (i in seq_len(nrow(M))) {
  b <- np$load(M$npz[i]); R_aug <- b[["R_aug"]]; z_aug <- as.numeric(b[["z_aug"]])
  t0 <- Sys.time()
  s <- tryCatch(estimate_s_rss(z = z_aug, R = R_aug), error = function(e) NA_real_)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  writeLines(paste(M$phys[i], M$trait[i], M$arm[i], M$region[i], M$reason[i], M$p[i],
                   sprintf("%.4f", s), sprintf("%.1f", wall), sep = "\t"), con)
  flush(con)
}
close(con)
cat(sprintf("[srss_run] %d cells -> %s\n", nrow(M), out))
