#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: run_comprehensive_susie_batch.R manifest.tsv variables.tsv status.tsv")
manifest_path <- args[[1]]
variables_path <- args[[2]]
status_path <- args[[3]]
manifest <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
variables <- list()
statuses <- list()
for (i in seq_len(nrow(manifest))) {
  pair_id <- manifest$pair_id[[i]]
  set.seed(as.integer(manifest$seed[[i]]))
  result <- tryCatch({
    dat <- read.delim(manifest$input_path[[i]], check.names = FALSE)
    y <- dat[["phenotype"]]
    X <- as.matrix(dat[, setdiff(names(dat), "phenotype"), drop = FALSE])
    if (ncol(X) < 1) stop("no predictors")
    fit <- susieR::susie(X, y, L = min(10L, ncol(X)), coverage = 0.95,
                         standardize = TRUE, intercept = FALSE, max_iter = 500,
                         estimate_residual_variance = TRUE)
    out <- data.frame(pair_id = pair_id, variable = colnames(X), pip = fit$pip,
                      credible_sets = "", stringsAsFactors = FALSE)
    if (length(fit$sets$cs) > 0) {
      for (j in seq_along(fit$sets$cs)) {
        idx <- fit$sets$cs[[j]]
        label <- paste0("CS", j)
        out$credible_sets[idx] <- ifelse(out$credible_sets[idx] == "", label,
                                         paste(out$credible_sets[idx], label, sep = ","))
      }
    }
    variables[[length(variables) + 1L]] <- out
    data.frame(pair_id = pair_id, status = ifelse(isTRUE(fit$converged), "ok", "not_converged"),
               converged = isTRUE(fit$converged), n_predictors = ncol(X),
               n_iterations = fit$niter, message = "", stringsAsFactors = FALSE)
  }, error = function(error) {
    data.frame(pair_id = pair_id, status = "failed", converged = FALSE,
               n_predictors = NA_integer_, n_iterations = NA_integer_,
               message = conditionMessage(error), stringsAsFactors = FALSE)
  })
  statuses[[length(statuses) + 1L]] <- result
}
variable_frame <- if (length(variables)) do.call(rbind, variables) else data.frame(
  pair_id = character(), variable = character(), pip = numeric(), credible_sets = character())
status_frame <- do.call(rbind, statuses)
write.table(variable_frame, variables_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(status_frame, status_path, sep = "\t", row.names = FALSE, quote = FALSE)
