# ============================================================================
# 08_Spatial_metabolic.R
# GSE181300 (8 Visium HNSCC): 标签转移 + 代谢共识基因空间 HubScore + GL/OM 空间评分
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

suppressPackageStartupMessages({library(Seurat); library(tidyverse); library(patchwork)})
set.seed(20260812)
ws <- root
ann_root <- hnscc_root
sp_base <- file.path(ann_root, "data", "spatial_samples")
dir.create(file.path(ws, "results", "Spatial"), showWarnings = FALSE, recursive = TRUE)
figdir <- file.path(ws, "figures")
samples <- list.files(sp_base); cat("spatial samples:", length(samples), "\n")

cons <- read.csv(file.path(ws, "results", "ML", "consensus_genes.csv"))
hub30 <- head(cons$gene, 30)
w <- setNames(head(cons$consensus_score, 30), hub30)
gs <- readRDS(file.path(ws, "results", "03_metabolic_genesets.rds"))

ref <- readRDS(file.path(ann_root, "results", "01_seurat_annotated.rds"))
ref <- DietSeurat(ref, assays = "RNA", layers = "data")
ref$cell_type <- factor(ref$cell_type)
ref <- FindVariableFeatures(ref, nfeatures = 3000, verbose = FALSE)
ref <- ScaleData(ref, verbose = FALSE)
ref <- RunPCA(ref, npcs = 30, verbose = FALSE)

all_spot <- list()
for (sname in samples) {
  cat("\n=== ", sname, " ===\n")
  s <- Load10X_Spatial(file.path(sp_base, sname),
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "RNA", image.name = "tissue_hires_image.png")
  s <- NormalizeData(s, verbose = FALSE)
  s <- FindVariableFeatures(s, nfeatures = 3000, verbose = FALSE)
  s <- ScaleData(s, verbose = FALSE)
  s <- RunPCA(s, npcs = 30, verbose = FALSE)

  anchors <- tryCatch(FindTransferAnchors(
    reference = ref, query = s, normalization.method = "LogNormalize",
    dims = 1:30, reduction = "pcaproject", verbose = FALSE),
    error = function(e) NULL)
  if (!is.null(anchors)) {
    td <- TransferData(anchors, refdata = ref$cell_type, dims = 1:30,
                       weight.reduction = "pcaproject", verbose = FALSE)
    s <- AddMetaData(s, td$predicted.id, "pred_celltype")
    s <- AddMetaData(s, td$prediction.score.max, "pred_score")
  }
  # 代谢 HubScore
  g_here <- intersect(hub30, rownames(s))
  expr <- as.matrix(s[["RNA"]]$data[g_here, , drop = FALSE])
  ww <- w[g_here]
  s$HubScore <- as.numeric(scale((t(expr) %*% ww) / sum(ww)))
  # GL/OM 空间代谢评分
  s <- AddModuleScore(s, features = list(intersect(gs$GL, rownames(s))),
                      name = "GLspatial", verbose = FALSE)
  s <- AddModuleScore(s, features = list(intersect(gs$OM, rownames(s))),
                      name = "OMspatial", verbose = FALSE)

  cc <- GetTissueCoordinates(s)
  df <- data.frame(sample = sname, spot = rownames(cc), x = cc$x, y = cc$y,
                   HubScore = s$HubScore, GLspatial = s$GLspatial1,
                   OMspatial = s$OMspatial1, pred_celltype = s$pred_celltype,
                   pred_score = s$pred_score)
  all_spot[[sname]] <- df
  write.csv(df, file.path(ws, "results", "Spatial", paste0(sname, "_spots.csv")),
            row.names = FALSE)

  p1 <- SpatialFeaturePlot(s, features = "HubScore", image.alpha = 0.6) +
    scale_fill_gradient2(low = "#123B7A", mid = "white", high = "#8A0612") +
    labs(title = paste(sname, "metabolic HubScore"))
  p_gl <- SpatialFeaturePlot(s, features = "GLspatial1", image.alpha = 0.6) +
    scale_fill_gradient(low = "white", high = "#8A0612") + labs(title="GL spatial")
  p_om <- SpatialFeaturePlot(s, features = "OMspatial1", image.alpha = 0.6) +
    scale_fill_gradient(low = "white", high = "#123B7A") + labs(title="OM spatial")
  pdf(file.path(figdir, paste0("08_", sname, "_spatial.pdf")), width = 13, height = 10)
  print(wrap_plots(list(p1, p_gl, p_om), ncol = 2))
  invisible(dev.off())
  rm(s, anchors); gc()
}
spot_all <- bind_rows(all_spot)
write.csv(spot_all, file.path(ws, "results", "Spatial", "all_spots_hubscore.csv"),
          row.names = FALSE)
cat("\ntotal spots:", nrow(spot_all), "\n")
summ <- spot_all %>% group_by(sample) %>%
  summarise(n_spots = n(), mean_Hub = mean(HubScore),
            GL_OM_cor = cor(GLspatial, OMspatial), .groups = "drop")
write.csv(summ, file.path(ws, "results", "Spatial", "sample_summary.csv"),
          row.names = FALSE)
print(summ)
cat("\nGL vs OM spatial Spearman overall:",
    cor(spot_all$GLspatial, spot_all$OMspatial, method = "spearman"), "\n")
p4 <- ggplot(spot_all, aes(sample, HubScore, fill = sample)) +
  geom_violin(trim = TRUE) + geom_boxplot(width = 0.12, outlier.shape = NA) +
  coord_flip() + theme_bw() + theme(legend.position = "none") +
  labs(x = NULL, y = "Spatial metabolic HubScore")
ggsave(file.path(figdir, "08_hubscore_violin.pdf"), p4, width = 8, height = 6)
ggsave(file.path(figdir, "08_hubscore_violin.png"), p4, width = 8, height = 6, dpi=300)
saveRDS(spot_all, file.path(ws, "results", "Spatial", "spatial_hubscore.rds"))
cat("\nDONE\n")
