# ============================================================================
# 01_QC_cluster.R
# GSE103322 HNSCC: QC -> HVG -> Scale -> PCA -> Harmony -> 聚类 -> UMAP
# 数据已是 log2(CP100K+1), 无 raw count, 无 MT 基因 -> 不做 NormalizeData, 不用 MT%
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
  library(Seurat); library(harmony); library(tidyverse); library(patchwork); library(ggsci)
})
set.seed(20260811)
root <- hnscc_root
fig <- file.path(root, "figures", "01_QC"); dir.create(fig, recursive = TRUE, showWarnings = FALSE)
res <- file.path(root, "results")

obj <- readRDS(file.path(res, "00_GSE103322_seurat_raw.rds"))
n0 <- ncol(obj)

# ---- 1. QC 门 (宽松; 数据已被原始质控) --------------------------------------
obj$QC_pass <- obj$nFeature_RNA >= 1500 & obj$nFeature_RNA <= 14000 &
               obj$nCount_RNA >= 3000 & obj$nCount_RNA <= 25000
cat("QC removed:", sum(!obj$QC_pass), "of", n0, "\n")
obj <- subset(obj, subset = QC_pass)

# ---- 2. 高变基因 (vst; counts 层即 log2 归一化值) --------------------------
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
top20 <- head(VariableFeatures(obj), 20)

# ---- 3. 缩放 + PCA ----------------------------------------------------------
obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, verbose = FALSE)

# ---- 4. Harmony 患者批次整合 ------------------------------------------------
obj <- RunHarmony(obj, group.by.vars = "patient", reduction = "pca",
                  reduction.save = "harmony", verbose = FALSE)

# ---- 5. 聚类 (多分辨率) + UMAP ----------------------------------------------
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
for (r in c(0.1, 0.2, 0.3, 0.5, 0.8)) {
  obj <- FindClusters(obj, resolution = r, verbose = FALSE, random.seed = 20260811)
}
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, seed.use = 2026811, verbose = FALSE)

# 主分辨率定为 0.3
obj$seurat_clusters <- obj$RNA_snn_res.0.3
Idents(obj) <- "RNA_snn_res.0.3"

saveRDS(obj, file.path(res, "01_seurat_clustered.rds"))
cat("clustered object saved:", ncol(obj), "cells,", length(levels(obj$RNA_snn_res.0.3)), "clusters @0.3\n")

# ---- 6. QC/UMAP 图 ----------------------------------------------------------
qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.ribo")
p1 <- VlnPlot(obj, features = qc_feats, group.by = "patient", pt.size = 0, ncol = 3) &
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
ggsave(file.path(fig, "01_qc_violin_by_patient.pdf"), p1, width = 13, height = 4.5)

p2 <- DimPlot(obj, reduction = "umap", group.by = "RNA_snn_res.0.3", label = TRUE) +
      ggtitle("Unsupervised clusters (res 0.3)")
p3 <- DimPlot(obj, reduction = "umap", group.by = "patient") + ggtitle("Patient")
p4 <- DimPlot(obj, reduction = "umap", group.by = "cancer_cell",
              cols = c("0" = "#4DBBD5", "1" = "#E64B35")) + ggtitle("Original cancer label")
ggsave(file.path(fig, "02_umap_cluster_patient_cancer.pdf"),
       p2 / (p3 | p4), width = 12, height = 11)

# PCA / Harmony 方差
el <- ElbowPlot(obj, ndims = 50, reduction = "pca")
ggsave(file.path(fig, "03_pca_elbow.pdf"), el, width = 5, height = 4)
cat("DONE\n")
