# ============================================================================
# 04c_visualize_pairing.R
# 转移配对可视化: 患者级配对连线(GL/OM) + 象限比例 + 关键基因配对
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
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
})
ws <- root
z <- readRDS(file.path(ws, "results", "03_malignant_grouped.rds"))
mal <- z$mal
paired_patients <- c("HN25","HN26","HN28","HNSCC20","HNSCC25",
                     "HNSCC26","HNSCC28","HNSCC5")
mp <- subset(mal, patient %in% paired_patients)
figdir <- file.path(ws, "figures")
theme_set(theme_classic(base_size = 12))
sitecol <- c("Primary" = "#123B7A", "Lymph node" = "#8A0612")

# ---- 患者级配对连线 --------------------------------------------------------
psum <- mp@meta.data %>% group_by(patient, lymph_node) %>%
  summarise(GL = mean(GL_score), OM = mean(OM_score), .groups = "drop") %>%
  mutate(site = ifelse(lymph_node == 1, "Lymph node", "Primary"),
         site = factor(site, levels = c("Primary","Lymph node")))
make_pair <- function(metric, ylab, pval) {
  ggplot(psum, aes(site, .data[[metric]], group = patient)) +
    geom_line(color = "grey50") +
    geom_point(aes(color = site), size = 2.5) +
    scale_color_manual(values = sitecol) +
    annotate("text", x = 1.5, y = Inf, label = paste0("p = ", pval),
             vjust = 2, size = 3.5) +
    labs(x = NULL, y = ylab, title = ylab) +
    theme(legend.position = "none")
}
p_gl <- make_pair("GL", "Glycolysis/Lactate score", "0.234")
p_om <- make_pair("OM", "OXPHOS score", "0.834")

# ---- 象限比例 --------------------------------------------------------------
qprop <- mp@meta.data %>% group_by(lymph_node, quadrant) %>%
  summarise(cells = n(), .groups = "drop") %>% group_by(lymph_node) %>%
  mutate(prop = cells / sum(cells),
         site = ifelse(lymph_node == 1, "Lymph node", "Primary"))
qcol <- c("Glycolytic_GLh_OMl" = "#8A0612", "Hybrid_GLh_OMh" = "#D8B96A",
          "Oxidative_GLl_OMh" = "#123B7A", "Quiescent_GLl_OMl" = "grey60")
p_q <- ggplot(qprop, aes(site, prop, fill = quadrant)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = qcol, labels = c("Glycolytic","Hybrid","Oxidative","Quiescent"),
                    name = "Quadrant") +
  scale_x_discrete(limits = c("Primary","Lymph node")) +
  labs(x = NULL, y = "Cell proportion", title = "Quadrant composition")

# ---- 关键基因配对 ----------------------------------------------------------
g <- c("GAPDH","LDHA","SLC16A3")
expr <- GetAssayData(mp, layer = "data")[g, , drop = FALSE]
gd <- mp@meta.data[, c("patient","lymph_node")]
glong <- bind_rows(lapply(g, function(gg) {
  data.frame(gene = gg, patient = gd$patient, lymph_node = gd$lymph_node,
             value = expr[gg, ]) %>% group_by(gene, patient, lymph_node) %>%
    summarise(v = mean(value), .groups = "drop")
})) %>% mutate(site = factor(ifelse(lymph_node == 1, "Lymph node", "Primary"),
                            levels = c("Primary","Lymph node")))
pvals <- c(GAPDH = "0.059", LDHA = "0.294", SLC16A3 = "0.294")
p_gene <- ggplot(glong, aes(site, v, group = patient)) +
  geom_line(color = "grey50") + geom_point(aes(color = site), size = 1.8) +
  scale_color_manual(values = sitecol) +
  facet_wrap(~ gene, scales = "free_y") +
  geom_text(data = data.frame(gene = g, lab = paste0("p=", pvals[g])),
            aes(x = 1.5, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.5, size = 3) +
  labs(x = NULL, y = "log2 expression", title = "Key glycolytic/lactate genes") +
  theme(legend.position = "none")

combined <- (p_gl | p_om) / p_q / p_gene +
  plot_annotation(title = "Primary vs lymph-node metastatic malignant cells (8 paired patients)")
ggsave(file.path(figdir, "04_transfer_pairing.pdf"), combined, width = 11, height = 14)
ggsave(file.path(figdir, "04_transfer_pairing.png"), combined, width = 11, height = 14, dpi = 300)
cat("已保存配对图\n")
