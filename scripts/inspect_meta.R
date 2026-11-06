# ---- Auto-locate and source project configuration (no hard-coded paths) --
..cfg_file <- (function() {
  d <- normalizePath(getwd(), mustWork = FALSE)
  for (i in seq_len(12)) {
    cand <- file.path(d, "config.R")
    if (file.exists(cand)) return(cand)
    up <- dirname(d)
    if (identical(up, d)) break
    d <- up
  }
  stop("config.R not found; please run from within the repository.")
})()
source(..cfg_file)

suppressPackageStartupMessages({library(Seurat); library(dplyr)})
f <- file.path(hnscc_root, "results", "01_seurat_annotated.rds")
obj <- readRDS(f)
cat("维度:", nrow(obj), "genes x", ncol(obj), "cells\n")
cat("=== meta列名 ===\n"); print(colnames(obj@meta.data))
cat("=== reductions ===\n"); print(Reductions(obj))
for (cc in colnames(obj@meta.data)) {
  v <- obj@meta.data[[cc]]
  if (is.character(v) || is.factor(v) || length(unique(v)) <= 12) {
    cat("\n--- ", cc, " (n=", length(unique(v)), ") ---\n", sep="")
    print(table(v, useNA="ifany"))
  } else {
    cat("\n--- ", cc, " 连续/高基数, 范围: ", min(v, na.rm=TRUE), "~", max(v, na.rm=TRUE), "\n", sep="")
  }
}
