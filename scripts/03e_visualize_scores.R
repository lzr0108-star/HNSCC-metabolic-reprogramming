# ============================================================================
# 03e_visualize_scores.R
# 代谢评分与分组可视化: UMAP(GL/OM/象限) + GL-OM 双轴散点
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

suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(patchwork); library(dplyr)})
ws <- root
z <- readRDS(file.path(ws, "results", "03_malignant_grouped.rds"))
mal <- z$mal
um <- as.data.frame(Embeddings(mal, "umap")); colnames(um) <- c("UMAP1","UMAP2")
um$GL <- mal$GL_score; um$OM <- mal$OM_score; um$quadrant <- mal$quadrant
um$subtype <- mal$cell_subtype
figdir <- file.path(ws, "figures")

col_hi <- "#8A0612"; col_lo <- "#123B7A"; col_mid <- "#D8B96A"
theme_set(theme_classic(base_size = 12))

p_gl <- ggplot(um, aes(UMAP1, UMAP2, color = GL)) +
  geom_point(size = 0.7, alpha = 0.85) +
  scale_color_gradient2(low = col_lo, mid = "white", high = col_hi, midpoint = 0) +
  labs(title = "Glycolysis/Lactate score", color = "GL")
p_om <- ggplot(um, aes(UMAP1, UMAP2, color = OM)) +
  geom_point(size = 0.7, alpha = 0.85) +
  scale_color_gradient2(low = col_lo, mid = "white", high = col_hi, midpoint = 0) +
  labs(title = "OXPHOS score", color = "OM")
qcol <- c("Glycolytic_GLh_OMl" = "#8A0612", "Hybrid_GLh_OMh" = "#D8B96A",
          "Oxidative_GLl_OMh" = "#123B7A", "Quiescent_GLl_OMl" = "grey60")
p_q <- ggplot(um, aes(UMAP1, UMAP2, color = quadrant)) +
  geom_point(size = 0.7, alpha = 0.85) +
  scale_color_manual(values = qcol, name = "Metabolic quadrant") +
  labs(title = "Metabolic quadrant") +
  theme(legend.position = "right")

p_scatter <- ggplot(um, aes(GL, OM, color = quadrant)) +
  geom_point(size = 0.7, alpha = 0.8) +
  geom_vline(xintercept = median(um$GL), linetype = 2, color = "grey40") +
  geom_hline(yintercept = median(um$OM), linetype = 2, color = "grey40") +
  scale_color_manual(values = qcol) +
  labs(x = "Glycolysis/Lactate score", y = "OXPHOS score",
       title = "Dual-axis metabolic landscape") +
  theme(legend.position = "none")

combined <- (p_gl | p_om) / (p_q | p_scatter) +
  plot_annotation(title = "HNSCC malignant cell metabolic heterogeneity (n=2500)")
ggsave(file.path(figdir, "03_metabolic_umap_quadrant.pdf"), combined,
       width = 11, height = 10)
ggsave(file.path(figdir, "03_metabolic_umap_quadrant.png"), combined,
       width = 11, height = 10, dpi = 300)

# 象限计数条形
p_count <- ggplot(um, aes(quadrant, fill = quadrant)) +
  geom_bar() + geom_text(stat = "count", aes(label = after_stat(count)),
                         vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = qcol) +
  scale_x_discrete(labels = c("Glycolytic","Hybrid","Oxidative","Quiescent")) +
  labs(x = NULL, y = "Malignant cells", title = "Quadrant cell counts") +
  theme(legend.position = "none")
ggsave(file.path(figdir, "03_quadrant_counts.pdf"), p_count, width = 6, height = 4.5)
ggsave(file.path(figdir, "03_quadrant_counts.png"), p_count, width = 6, height = 4.5, dpi=300)
cat("已保存评分与分组图\n")
