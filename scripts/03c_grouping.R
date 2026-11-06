# ============================================================================
# 03c_grouping.R
# 方案A: GL/OM 双轴中位数四象限; 方案B: GL 单轴 Q25/Q75 三分组
# 与恶性亚型、lymph_node 交叉; 输出分组表与更新对象
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

suppressPackageStartupMessages({library(Seurat); library(dplyr)})
ws <- root
z <- readRDS(file.path(ws, "results", "03_malignant_metabolic_scores.rds"))
mal <- z$mal; score_df <- z$score_df

# ---- 方案 A: 双轴中位数四象限 ----------------------------------------------
gl_med <- median(score_df$GL_score); om_med <- median(score_df$OM_score)
score_df$GL_high <- score_df$GL_score >= gl_med
score_df$OM_high <- score_df$OM_score >= om_med
score_df$quadrant <- ifelse(score_df$GL_high & score_df$OM_high, "Hybrid_GLh_OMh",
                    ifelse(score_df$GL_high & !score_df$OM_high, "Glycolytic_GLh_OMl",
                    ifelse(!score_df$GL_high & score_df$OM_high, "Oxidative_GLl_OMh",
                    "Quiescent_GLl_OMl")))
cat("=== 方案A 四象限分布 ===\n"); print(table(score_df$quadrant))

# ---- 方案 B: GL 单轴四分位 --------------------------------------------------
q25 <- quantile(score_df$GL_score, 0.25); q75 <- quantile(score_df$GL_score, 0.75)
score_df$GL_group <- ifelse(score_df$GL_score >= q75, "HighGL",
                     ifelse(score_df$GL_score <= q25, "LowGL", "DTGL"))
score_df$GL_group <- factor(score_df$GL_group, levels = c("LowGL","DTGL","HighGL"))
cat("\n=== 方案B GL三分组 ===\n"); print(table(score_df$GL_group))

# OM 单轴四分位(配套, 便于后续)
o25 <- quantile(score_df$OM_score, 0.25); o75 <- quantile(score_df$OM_score, 0.75)
score_df$OM_group <- ifelse(score_df$OM_score >= o75, "HighOM",
                     ifelse(score_df$OM_score <= o25, "LowOM", "DTOM"))
score_df$OM_group <- factor(score_df$OM_group, levels = c("LowOM","DTOM","HighOM"))

# ---- 与恶性亚型交叉 --------------------------------------------------------
cat("\n=== 四象限 x 恶性亚型 ===\n")
print(table(score_df$quadrant, score_df$cell_subtype))
cat("\n=== 四象限 lymph_node 构成 ===\n")
print(prop.table(table(score_df$quadrant, score_df$lymph_node), margin = 1))
cat("\n=== GL组 lymph_node 构成 ===\n")
print(prop.table(table(score_df$GL_group, score_df$lymph_node), margin = 1))

# ---- 写回 Seurat 元数据 ----------------------------------------------------
ix <- match(rownames(mal@meta.data), score_df$cell_id)
mal$quadrant <- score_df$quadrant[ix]
mal$GL_group <- score_df$GL_group[ix]
mal$OM_group <- score_df$OM_group[ix]
mal$GL_score <- score_df$GL_score[ix]
mal$OM_score <- score_df$OM_score[ix]

write.csv(score_df, file.path(ws, "results", "03_metabolic_scores_grouped.csv"), row.names = FALSE)
saveRDS(list(mal = mal, score_df = score_df), file.path(ws, "results", "03_malignant_grouped.rds"))
cat("\n已保存分组结果\n")
cat("中位数切点 GL:", round(gl_med,3), "OM:", round(om_med,3),
    "| Q25/Q75 GL:", round(q25,3), round(q75,3), "\n")
