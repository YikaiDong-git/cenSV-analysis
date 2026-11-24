#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
input <- args[[1]]; output <- args[[2]]
set.seed(as.integer(args[[3]]))
dat <- read.delim(input, check.names = FALSE)
y <- dat[["phenotype"]]
X <- as.matrix(dat[, setdiff(names(dat), "phenotype"), drop = FALSE])
fit <- susieR::susie(X, y, L = min(10, ncol(X)), coverage = 0.95,
                     standardize = TRUE, intercept = TRUE, max_iter = 500)
out <- data.frame(variable = colnames(X), pip = fit$pip)
out$credible_sets <- ""
if (length(fit$sets$cs) > 0) {
  for (i in seq_along(fit$sets$cs)) {
    idx <- fit$sets$cs[[i]]
    out$credible_sets[idx] <- ifelse(out$credible_sets[idx] == "", paste0("CS", i),
                                     paste(out$credible_sets[idx], paste0("CS", i), sep = ","))
  }
}
write.table(out, output, sep = "\t", row.names = FALSE, quote = FALSE)
