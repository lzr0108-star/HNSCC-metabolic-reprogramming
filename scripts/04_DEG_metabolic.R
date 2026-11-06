# ============================================================================
# 04_DEG_metabolic.R
# 高糖酵解(HighGL) vs 低糖酵解(LowGL) 恶性细胞差异表达
# 输入: 03_malignant_grouped.rds; Wilcox, BH FDR<0.05 且 |log2FC|>=0.5
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

suppressPackageStartupMessages({library(Seurat); library(tidyverse)})
set.seed(20260811)
ws <- root
dir.create(file.path(ws, "results", "DEG"), showWarnings = FALSE, recursive = TRUE)
figdir <- file.path(ws, "figures"); dir.create(figdir, showWarnings = FALSE)

z <- readRDS(file.path(ws, "results", "03_malignant_grouped.rds"))
mal <- z$mal
Idents(mal) <- "GL_group"

deg <- FindMarkers(mal, ident.1 = "HighGL", ident.2 = "LowGL",
                   test.use = "wilcox", logfc.threshold = 0, min.pct = 0.1,
                   verbose = FALSE)
deg$gene <- rownames(deg)
deg <- deg[, c("gene","avg_log2FC","pct.1","pct.2","p_val","p_val_adj")]
deg <- deg[order(deg$p_val_adj), ]
write.csv(deg, file.path(ws, "results", "DEG", "DEG_HighGL_vs_LowGL.csv"),
          row.names = FALSE)

sig <- deg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) >= 0.5)
up <- sig %>% filter(avg_log2FC > 0); dn <- sig %>% filter(avg_log2FC < 0)
write.csv(up, file.path(ws, "results", "DEG", "DEG_HighGL_up.csv"), row.names = FALSE)
write.csv(dn, file.path(ws, "results", "DEG", "DEG_HighGL_down.csv"), row.names = FALSE)
cat("total tested:", nrow(deg), "| significant:", nrow(sig),
    "| up:", nrow(up), "| down:", nrow(dn), "\n")
cat("Top 20 up:\n"); print(head(up$gene, 20))
cat("Top 20 down:\n"); print(head(dn$gene, 20))

# ---- 代谢基因方向自洽性 ----------------------------------------------------
gs <- readRDS(file.path(ws, "results", "03_metabolic_genesets.rds"))
check <- deg %>% filter(gene %in% c(gs$GL, gs$OM)) %>%
  mutate(axis = ifelse(gene %in% gs$GL, "GL_gene", "OM_gene"))
cat("\n=== 评分基因在 HighGL 对比中的 logFC 方向 ===\n")
print(check %>% group_by(axis) %>% summarise(
  n = n(), median_logFC = median(avg_log2FC),
  n_pos = sum(avg_log2FC > 0), n_neg = sum(avg_log2FC < 0)))

# ---- 火山图 ----------------------------------------------------------------
deg$class <- ifelse(deg$p_val_adj < 0.05 & deg$avg_log2FC >= 0.5, "Up (HighGL)",
              ifelse(deg$p_val_adj < 0.05 & deg$avg_log2FC <= -0.5, "Down (LowGL)", "NS"))
p1 <- ggplot(deg, aes(avg_log2FC, -log10(p_val_adj), colour = class)) +
  geom_point(size = 0.7, alpha = 0.7) +
  scale_colour_manual(values = c("Up (HighGL)" = "#8A0612",
                                 "Down (LowGL)" = "#123B7A", "NS" = "grey80")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey50") +
  theme_bw() + labs(x = "avg log2FC (HighGL vs LowGL)", y = "-log10(adj p)", colour = NULL)
ggsave(file.path(figdir, "04_DEG_volcano.pdf"), p1, width = 6.5, height = 5)
ggsave(file.path(figdir, "04_DEG_volcano.png"), p1, width = 6.5, height = 5, dpi = 300)

# ---- Top DEG 热图 ----------------------------------------------------------
topg <- head(c(up$gene, dn$gene), 40)
cells <- colnames(mal)[mal$GL_group %in% c("HighGL","LowGL")]
p2 <- DoHeatmap(subset(mal, cells = cells), features = topg, group.by = "GL_group",
                group.colors = c("HighGL" = "#8A0612", "LowGL" = "#123B7A")) +
  scale_fill_gradient2(low = "#123B7A", mid = "white", high = "#8A0612")
ggsave(file.path(figdir, "04_DEG_heatmap.pdf"), p2, width = 9, height = 7)
ggsave(file.path(figdir, "04_DEG_heatmap.png"), p2, width = 9, height = 7, dpi = 300)
cat("\nDONE\n")
