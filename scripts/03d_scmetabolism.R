# ============================================================================
# 03d_scmetabolism.R
# 用 scMetabolism 内置 KEGG gmt(85 通路) + 新版 GSVA ssgseaParam 实现等价全景评分
# 与自建 GL/OM 综合评分做一致性验证
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
  library(Seurat); library(dplyr); library(GSEABase); library(GSVA)
})
set.seed(20260811)
ws <- root
z <- readRDS(file.path(ws, "results", "03_malignant_grouped.rds"))
mal <- z$mal; score_df <- z$score_df
expr <- as.matrix(GetAssayData(mal, layer = "data"))

gmt <- getGmt(system.file("data", "KEGG_metabolism_nc.gmt", package = "scMetabolism"))
# 过滤为对象内基因(显式保留通路名)
geneSets <- setNames(lapply(gmt, function(g) intersect(geneIds(g), rownames(expr))),
                     names(gmt))
geneSets <- geneSets[sapply(geneSets, length) >= 5]
cat("有效通路数:", length(geneSets), "\n"); print(head(names(geneSets), 4))

param <- ssgseaParam(expr, geneSets = geneSets, normalize = TRUE, verbose = FALSE)
kegg_scores <- gsva(param, verbose = FALSE)
cat("全景评分矩阵:", nrow(kegg_scores), "x", ncol(kegg_scores), "\n")

key <- rownames(kegg_scores)[grep("Glycolysis|TCA cycle|Pyruvate metabolism|Oxidative phosphorylation",
                                   rownames(kegg_scores))]
cat("匹配到关键通路:", length(key), "\n"); print(key)
cat("\n=== 关键通路与自建综合评分 Spearman ===\n")
for (k in key) {
  v <- kegg_scores[k, ]
  cat(sprintf("%-34s vs GL: %.3f  vs OM: %.3f\n", k,
              cor(v, score_df$GL_score, method = "spearman"),
              cor(v, score_df$OM_score, method = "spearman")))
}
# 关键通路间相关
cat("\n=== 关键通路间 Spearman ===\n")
print(round(cor(t(kegg_scores[key, ]), method = "spearman"), 3))

saveRDS(kegg_scores, file.path(ws, "results", "03_scMetabolism_KEGG_scores.rds"))
write.csv(as.data.frame(t(kegg_scores)),
          file.path(ws, "results", "03_scMetabolism_KEGG_scores.csv"))
cat("\n已保存全景代谢评分\n")
