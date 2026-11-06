# ============================================================================
# 03b_multi_scoring.R
# 恶性细胞 GL/OM 双轴多算法评分: AddModuleScore / AUCell / ssgsea(GSVA) / singscore
# 输出: results/03_malignant_metabolic_scores.rds ; 评分CSV
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
  library(Seurat); library(dplyr); library(Matrix)
  library(AUCell); library(GSVA); library(singscore)
})
set.seed(20260811)
ws <- root
gs <- readRDS(file.path(ws, "results", "03_metabolic_genesets.rds"))
GL <- gs$GL; OM <- gs$OM
sets <- list(GL = GL, OM = OM)

obj <- readRDS(file.path(hnscc_root, "results", "01_seurat_annotated.rds"))
mal <- subset(obj, cell_type == "Malignant")
rm(obj); invisible(gc())
cat("恶性细胞数:", ncol(mal), "\n")
expr <- as.matrix(GetAssayData(mal, layer = "data"))
cat("表达矩阵:", nrow(expr), "x", ncol(expr), "\n")
score_df <- data.frame(cell_id = colnames(mal), stringsAsFactors = FALSE)

# ---- 1. AddModuleScore -----------------------------------------------------
mal <- AddModuleScore(mal, features = list(GL = GL), name = "GL_ams")
mal <- AddModuleScore(mal, features = list(OM = OM), name = "OM_ams")
score_df$GL_ams <- mal$GL_ams1; score_df$OM_ams <- mal$OM_ams1
cat("AddModuleScore 完成\n")

# ---- 2. AUCell --------------------------------------------------------------
rankings <- AUCell_buildRankings(expr, plotStats = FALSE, verbose = FALSE)
auc <- AUCell_calcAUC(sets, rankings, verbose = FALSE)
aucmat <- getAUC(auc)
score_df$GL_aucell <- aucmat["GL", ]; score_df$OM_aucell <- aucmat["OM", ]
cat("AUCell 完成\n")

# ---- 3. GSVA ssgsea ---------------------------------------------------------
param <- ssgseaParam(expr, geneSets = sets, normalize = TRUE, verbose = FALSE)
ss <- gsva(param, verbose = FALSE)
score_df$GL_ssgsea <- ss["GL", ]; score_df$OM_ssgsea <- ss["OM", ]
cat("ssgsea 完成\n")

# ---- 4. singscore -----------------------------------------------------------
rankData <- rankGenes(expr)
gl_s <- simpleScore(rankData, upSet = GL)
om_s <- simpleScore(rankData, upSet = OM)
score_df$GL_singscore <- gl_s$TotalScore
score_df$OM_singscore <- om_s$TotalScore
cat("singscore 完成\n")

# ---- 跨算法 Spearman 一致性 ------------------------------------------------
algs <- c("ams","aucell","ssgsea","singscore")
cat("\n=== GL 跨算法 Spearman ===\n")
gl_cols <- paste0("GL_", algs)
print(round(cor(score_df[gl_cols], method = "spearman"), 3))
cat("=== OM 跨算法 Spearman ===\n")
om_cols <- paste0("OM_", algs)
print(round(cor(score_df[om_cols], method = "spearman"), 3))

# ---- 综合评分: 4 算法 z-score 均值 -----------------------------------------
z <- function(x) as.numeric(scale(x))
score_df$GL_score <- rowMeans(sapply(score_df[gl_cols], z))
score_df$OM_score <- rowMeans(sapply(score_df[om_cols], z))
# 纳入关键元数据
md <- mal@meta.data
keep_cols <- c("patient","lymph_node","cancer_cell","is_combo","cell_subtype","seurat_clusters")
score_df <- cbind(score_df, md[match(score_df$cell_id, rownames(md)), keep_cols, drop=FALSE])

write.csv(score_df, file.path(ws, "results", "03_metabolic_scores.csv"), row.names = FALSE)
saveRDS(list(mal = mal, score_df = score_df), file.path(ws, "results", "03_malignant_metabolic_scores.rds"))
cat("\n已保存评分, 细胞数:", nrow(score_df), "\n")
