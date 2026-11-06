# ===========================================================================
# Project configuration
# Reproducible code for the HNSCC glycolysis/lactate (GL) vs mitochondrial
# OXPHOS (OM) metabolic reprogramming study.
# This file is located automatically and sourced by every analysis script,
# so no hard-coded absolute paths are needed.
# ===========================================================================

# ---- Locate the repository root (it contains 'scripts/' and this file) ---
.find_repo_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = FALSE)
  for (i in seq_len(12)) {
    if (dir.exists(file.path(d, "scripts")) &&
        file.exists(file.path(d, "config.R"))) {
      return(d)
    }
    up <- dirname(d)
    if (identical(up, d)) break
    d <- up
  }
  stop("Repository root not found. Please run from within the repository.")
}

# ---- Standard directory handles for this repository ----------------------
root    <- .find_repo_root()
res_dir <- file.path(root, "results")
fig_dir <- file.path(root, "figures")
tom_dir <- file.path(root, "TOM")

# ---- External source project ------------------------------------------------
# Prerequisite steps 00-02 (object build, QC/clustering/annotation, CopyKAT),
# the raw Visium spatial samples and the TCGA GDC inputs live in the original
# HNSCC project. Point to it either with the HNSCC_PROJECT_ROOT environment
# variable, or by placing that project as a sibling folder named "HNSCC复现".
hnscc_root <- Sys.getenv("HNSCC_PROJECT_ROOT", unset = "")
if (!nzchar(hnscc_root)) {
  hnscc_root <- file.path(dirname(root), "HNSCC复现")
}

# ---- Convenience helper: ensure a directory exists -----------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

cat("Configuration loaded. Repository root:", root, "\n")
