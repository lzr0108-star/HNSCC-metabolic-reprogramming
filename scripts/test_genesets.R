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

suppressPackageStartupMessages({library(Seurat); library(dplyr); library(msigdbr)})
obj <- readRDS(file.path(hnscc_root, "results", "01_seurat_annotated.rds"))
mal <- subset(obj, cell_type == "Malignant")
cat("恶性细胞数:", ncol(mal), "\n")
cat("=== 恶性细胞 lymph_node ===\n"); print(table(mal$lymph_node))
cat("=== 恶性细胞 patient x lymph_node ===\n")
tab <- table(mal$patient, mal$lymph_node)
print(tab)
paired <- rownames(tab)[tab[,1] > 0 & tab[,2] > 0]
cat("配对患者数:", length(paired), "\n"); print(paired)

cat("\n=== 测试 msigdbr ===\n")
hm <- msigdbr(species = "Homo sapiens", category = "H")
cat("Hallmark基因集数:", length(unique(hm$gs_name)), "\n")
gl <- hm %>% filter(gs_name == "HALLMARK_GLYCOLYSIS") %>% pull(gene_symbol)
ox <- hm %>% filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>% pull(gene_symbol)
cat("HALLMARK_GLYCOLYSIS 基因数:", length(gl), "\n")
cat("HALLMARK_OXPHOS 基因数:", length(ox), "\n")
cat("两集合重叠:", length(intersect(gl, ox)), "\n")
