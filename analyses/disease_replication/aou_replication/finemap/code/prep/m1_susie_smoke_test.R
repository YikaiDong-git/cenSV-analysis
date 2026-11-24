#!/usr/bin/env Rscript
.libPaths(c(path.expand("~/.local/lib/R"), .libPaths()))
suppressMessages(library(susieR))
cat("susieR version:", as.character(packageVersion("susieR")), "\n\n")

in_cs <- function(csl, v) length(csl) > 0 && any(vapply(csl, function(s) v %in% s, logical(1)))

check <- function(name, truth) {
  e <- new.env(); data(list = name, package = "susieR", envir = e); d <- get(name, envir = e)
  y <- d$Y[, 1]
  f_ind <- suppressWarnings(susie(d$X, y, L = 10))
  R <- cor(d$X)
  r <- as.vector(cor(d$X, y)); r <- pmin(pmax(r, -0.999999), 0.999999)
  z <- r * sqrt((nrow(d$X) - 2) / (1 - r^2))
  f_rss <- suppressWarnings(susie_rss(z, R, n = nrow(d$X), L = 10))
  cs_i <- if (is.null(f_ind$sets$cs)) list() else f_ind$sets$cs
  cs_r <- if (is.null(f_rss$sets$cs)) list() else f_rss$sets$cs
  cap_i <- all(vapply(truth, function(v) in_cs(cs_i, v), logical(1)))
  cap_r <- all(vapply(truth, function(v) in_cs(cs_r, v), logical(1)))
  pc <- cor(f_ind$pip, f_rss$pip)
  cat(sprintf("%-14s truth={%s} | individual: %d CS, captured=%s | rss: %d CS, captured=%s | PIPcor(ind,rss)=%.4f\n",
              name, paste(truth, collapse = ","), length(cs_i), cap_i, length(cs_r), cap_r, pc))
  c(cap_i = cap_i, cap_r = cap_r, pc = pc)
}

a <- check("N3finemapping", c(403, 653, 773))
b <- check("N2finemapping", c(337, 999))
pass <- all(c(a["cap_i"], a["cap_r"], b["cap_i"], b["cap_r"]) == 1) && a["pc"] > 0.99 && b["pc"] > 0.99
cat("\n", if (pass) "PASS - engine recovers both known answers; individual~rss PIP corr > 0.99" else "FAIL - investigate", "\n")
quit(status = if (pass) 0L else 1L)
