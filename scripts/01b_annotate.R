# ============================================================================
# 01b_annotate.R
# 基于 canonical marker + 原始注释的细胞类型注释 (两级: 大类 + 亚型)
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

suppressPackageStartupMessages({ library(Seurat); library(tidyverse); library(ggsci) })
set.seed(20260811)
root <- hnscc_root
fig <- file.path(root, "figures", "01_QC")
obj <- readRDS(file.path(root, "results", "01_seurat_clustered.rds"))
cl <- droplevels(as.factor(obj$RNA_snn_res.0.3))

# ---- cluster -> 亚型 / 大类 -------------------------------------------------
sub_map <- c(
  "0"="Malignant differentiated", "1"="CAF", "2"="CAF",
  "3"="Cytotoxic T/NK", "4"="Malignant cycling", "5"="Treg",
  "6"="Malignant ALDH/metabolic", "7"="Malignant pEMT", "8"="Endothelial",
  "9"="Malignant cycling", "10"="Plasma", "11"="Mast",
  "12"="Myocyte/myofibroblast", "13"="Macrophage", "14"="Dendritic",
  "15"="Endothelial", "16"="B cell"
)
major_map <- c(
  "0"="Malignant", "1"="CAF", "2"="CAF", "3"="T/NK", "4"="Malignant",
  "5"="Treg", "6"="Malignant", "7"="Malignant", "8"="Endothelial",
  "9"="Malignant", "10"="Plasma", "11"="Mast", "12"="Myocyte",
  "13"="Macrophage", "14"="Dendritic", "15"="Endothelial", "16"="B cell"
)
obj$cell_subtype <- factor(unname(sub_map[as.character(cl)]))
obj$cell_type    <- factor(unname(major_map[as.character(cl)]))
Idents(obj) <- "cell_type"

# ---- 与原始 cancer 标签一致性核对 ------------------------------------------
cat("=== major cell type x original cancer label ===\n")
print(table(obj$cell_type, obj$cancer_cell))
cat("\n=== cell type counts ===\n"); print(table(obj$cell_type))
mal <- obj$cell_type == "Malignant"
agree <- mean(mal == (obj$cancer_cell == 1))
cat("\nagreement with original cancer label:", round(agree, 4), "\n")

# ---- 注释 UMAP 图 -----------------------------------------------------------
pal <- c("Malignant"="#8A0612", "CAF"="#D8B96A", "T/NK"="#123B7A", "Treg"="#7E6148",
         "B cell"="#00A087", "Plasma"="#3C5488", "Endothelial"="#E64B35",
         "Macrophage"="#F39B7F", "Mast"="#B09C85", "Dendritic"="#8491B4", "Myocyte"="#91D1C2")
p1 <- DimPlot(obj, group.by = "cell_type", cols = pal, label = FALSE) + ggtitle("Cell type")
p2 <- DimPlot(obj, group.by = "cell_subtype", label = TRUE) + ggtitle("Subtype") +
      theme(legend.position = "none")
ggsave(file.path(fig, "05_umap_annotated.pdf"), p1 | p2, width = 15, height = 6.5)

# 按患者/淋巴结的组成
comp <- obj@meta.data %>% filter(!is_combo) %>%
  count(patient, cell_type) %>% group_by(patient) %>% mutate(prop = n/sum(n))
p3 <- ggplot(comp, aes(patient, prop, fill = cell_type)) +
  geom_col() + scale_fill_manual(values = pal) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text = element_text(size = 7)) +
  labs(x = "", y = "Proportion", fill = "Cell type")
ggsave(file.path(fig, "06_composition_by_patient.pdf"), p3, width = 10, height = 5)

saveRDS(obj, file.path(root, "results", "01_seurat_annotated.rds"))
cat("\nsaved annotated object. DONE\n")
