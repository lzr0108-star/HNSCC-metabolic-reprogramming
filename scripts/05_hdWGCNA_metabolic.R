# ============================================================================
# 05_hdWGCNA_metabolic.R
# 恶性细胞 hdWGCNA 共表达网络 -> 代谢(GL/OM)相关模块 + Hub 基因
# 网络参数复用: 3000 HVG 去 RPS/RPL/HB/MT; metacell k=20; power 8; signed; bicor
# 模块-性状改为 GL_score/OM_score/lymph_node/cell_subtype
# ============================================================================
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

suppressPackageStartupMessages({
  library(Seurat); library(hdWGCNA); library(patchwork); library(tidyverse)
})
set.seed(20260811)
ws <- root
dir.create(file.path(ws, "results", "hdWGCNA"), showWarnings = FALSE, recursive = TRUE)
figdir <- file.path(ws, "figures")

z <- readRDS(file.path(ws, "results", "03_malignant_grouped.rds"))
mal <- z$mal
# 确保有 VariableFeatures
if (length(VariableFeatures(mal)) < 500) {
  mal <- FindVariableFeatures(mal, nfeatures = 3000, verbose = FALSE)
}
feat <- VariableFeatures(mal)
feat <- feat[!grepl("^RP[SL]|^HB[AB]|^MT-", feat)]
cat("network genes:", length(feat), "\n")
mal <- SetupForWGCNA(mal, wgcna_name = "metabolic", features = feat)

# ---- Metacell (k=20, 按亚型) -----------------------------------------------
mal <- MetacellsByGroups(
  seurat_obj = mal, group.by = "cell_subtype", ident.group = "cell_subtype",
  k = 20, reduction = "harmony", dims = 1:30,
  max_shared = 10, min_cells = 20, verbose = FALSE)
mal <- NormalizeMetacells(mal); mal <- ScaleMetacells(mal)
cat("metacells:", ncol(GetMetacellObject(mal)), "\n")

# ---- soft power ------------------------------------------------------------
mal <- SetDatExpr(mal, group_name = NULL, use_metacells = TRUE)
mal <- TestSoftPowers(mal, networkType = "signed", corFnc = "bicor")
p_sp <- PlotSoftPowers(mal)
pdf(file.path(figdir, "05_soft_power.pdf"), width = 9, height = 6)
wrap_plots(p_sp, ncol = 2); invisible(dev.off())

# ---- 构建网络 (power 8) -----------------------------------------------------
mal <- ConstructNetwork(
  mal, soft_power = 8, networkType = "signed", TOMType = "signed",
  corType = "bicor", TOMDenom = "min",
  minModuleSize = 30, mergeCutHeight = 0.2, deepSplit = 4,
  pamStage = FALSE, maxBlockSize = 30000, overwrite_tom = TRUE)
mods <- GetModules(mal); mod_tab <- table(mods$module)
cat("=== module sizes ===\n"); print(mod_tab)
write.csv(as.data.frame(mod_tab), file.path(ws, "results", "hdWGCNA", "module_sizes.csv"),
          row.names = FALSE)

# ---- MEs + kME -------------------------------------------------------------
mal <- ModuleEigengenes(mal, harmonized = FALSE, verbose = FALSE)
mal <- ModuleConnectivity(mal, harmonized = FALSE, corFnc = "bicor",
                          reassign_modules = FALSE)
mods <- GetModules(mal)
write.csv(mods, file.path(ws, "results", "hdWGCNA", "gene_module_kME.csv"),
          row.names = FALSE)
pdf(file.path(figdir, "05_dendrogram.pdf"), width = 10, height = 5)
PlotDendrogram(mal, main = "Metabolic hdWGCNA dendrogram"); invisible(dev.off())

# ---- module UMAP + hubs -----------------------------------------------------
mal <- RunModuleUMAP(mal, n_hubs = 10, exclude_grey = TRUE, n_neighbors = 15,
                     min_dist = 0.3)
p_umap <- ModuleUMAPPlot(mal, point_size = 1.2, label_hubs = 3)
ggsave(file.path(figdir, "05_module_umap.pdf"), p_umap, width = 8, height = 6.5)
ggsave(file.path(figdir, "05_module_umap.png"), p_umap, width = 8, height = 6.5, dpi=300)

hub <- GetHubGenes(mal, n_hubs = 10)
write.csv(hub, file.path(ws, "results", "hdWGCNA", "hub_genes.csv"), row.names = FALSE)

# ---- 模块-性状相关性 (代谢性状) --------------------------------------------
mal <- ModuleTraitCorrelation(
  mal, traits = c("GL_score","OM_score","lymph_node","cell_subtype"),
  group.by = "cell_subtype", cor_method = "pearson")
mtc <- GetModuleTraitCorrelation(mal)
cor_mat <- mtc$cor$all_cells; p_mat <- mtc$pval$all_cells; fdr_mat <- mtc$fdr$all_cells
write.csv(cor_mat, file.path(ws, "results", "hdWGCNA", "module_trait_cor.csv"))
write.csv(p_mat, file.path(ws, "results", "hdWGCNA", "module_trait_p.csv"))
write.csv(fdr_mat, file.path(ws, "results", "hdWGCNA", "module_trait_fdr.csv"))
p_mt <- PlotModuleTraitCorrelation(mal)
ggsave(file.path(figdir, "05_module_trait_heatmap.pdf"), p_mt, width = 9, height = 6.5)
ggsave(file.path(figdir, "05_module_trait_heatmap.png"), p_mt, width = 9, height = 6.5, dpi=300)

cat("\n=== module correlation with GL score ===\n")
print(round(sort(cor_mat["GL_score", ], decreasing = TRUE), 3))
cat("\n=== module correlation with OM score ===\n")
print(round(sort(cor_mat["OM_score", ], decreasing = TRUE), 3))

saveRDS(mal, file.path(ws, "results", "05_metabolic_hdWGCNA.rds"))
cat("\nDONE\n")
