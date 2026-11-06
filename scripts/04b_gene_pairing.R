# ============================================================================
# 04b_gene_pairing.R
# 关键糖酵解/乳酸/OXPHOS 单基因 患者级 原发 vs 淋巴结 配对检验
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

suppressPackageStartupMessages({library(Seurat); library(dplyr); library(tidyr)})
ws <- root
z <- readRDS(file.path(ws, "results", "03_malignant_grouped.rds"))
mal <- z$mal
paired_patients <- c("HN25","HN26","HN28","HNSCC20","HNSCC25",
                     "HNSCC26","HNSCC28","HNSCC5")
mp <- subset(mal, patient %in% paired_patients)
genes <- c("LDHA","LDHB","SLC16A1","SLC16A3","SLC16A7","HK2","PKM","PFPKP",
           "ENO1","GAPDH","PGK1","ALDOA",
           "CS","IDH2","FH","MDH2","COX5A","NDUFA1","ATP5F1A","SDHA")
genes <- genes[genes %in% rownames(mp)]
expr <- GetAssayData(mp, layer = "data")[genes, , drop = FALSE]
md <- mp@meta.data[, c("patient","lymph_node")]

# 患者 x 部位 基因均值
long <- lapply(genes, function(g) {
  data.frame(gene = g, patient = md$patient, lymph_node = md$lymph_node,
             value = expr[g, ]) %>% group_by(gene, patient, lymph_node) %>%
    summarise(v = mean(value), .groups = "drop")
}) %>% bind_rows()
w <- long %>% pivot_wider(names_from = lymph_node, values_from = v,
                          names_prefix = "site")
res <- w %>% group_by(gene) %>% summarise(
  mean_PR = mean(site0), mean_LN = mean(site1),
  delta = mean(site1 - site0), n_up = sum(site1 > site0),
  cohen_d = mean(site1 - site0) / sd(site1 - site0),
  p_value = wilcox.test(site1, site0, paired = TRUE, exact = FALSE)$p.value,
  .groups = "drop") %>% arrange(p_value)
print(as.data.frame(res), row.names = FALSE)
write.csv(res, file.path(ws, "results", "04_gene_paired_tests.csv"), row.names = FALSE)
