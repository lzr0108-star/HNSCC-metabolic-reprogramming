# ============================================================================
# 04_transfer_pairing.R
# 原发(lymph_node=0) vs 淋巴结转移(lymph_node=1) 恶性细胞代谢配对分析
# 统计单位=患者(8 个配对), 配对 Wilcoxon; 避免把细胞当独立重复
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

suppressPackageStartupMessages({library(dplyr); library(tidyr)})
ws <- root
d <- read.csv(file.path(ws, "results", "03_metabolic_scores_grouped.csv"),
              stringsAsFactors = FALSE)
paired_patients <- c("HN25","HN26","HN28","HNSCC20","HNSCC25",
                     "HNSCC26","HNSCC28","HNSCC5")
dp <- d %>% filter(patient %in% paired_patients)
cat("配对分析细胞数:", nrow(dp), "| 患者数:", length(unique(dp$patient)), "\n")
cat("各患者原发/淋巴结细胞数:\n"); print(table(dp$patient, dp$lymph_node))

# ---- 患者 x 部位 均值汇总 --------------------------------------------------
agg <- dp %>% group_by(patient, lymph_node) %>% summarise(
  n = n(),
  GL = mean(GL_score), OM = mean(OM_score),
  GL_ams = mean(GL_ams), OM_ams = mean(OM_ams),
  GL_aucell = mean(GL_aucell), OM_aucell = mean(OM_aucell),
  GL_ssgsea = mean(GL_ssgsea), OM_ssgsea = mean(OM_ssgsea),
  GL_singscore = mean(GL_singscore), OM_singscore = mean(OM_singscore),
  .groups = "drop")
write.csv(agg, file.path(ws, "results", "04_patient_site_summary.csv"), row.names = FALSE)

# 宽格式: 每个患者一列原发、一列淋巴结
w <- agg %>% mutate(site = ifelse(lymph_node == 1, "LN", "PR")) %>%
  select(-lymph_node) %>%
  pivot_wider(names_from = site, values_from = c(n, GL, OM, GL_ams, OM_ams,
    GL_aucell, OM_aucell, GL_ssgsea, OM_ssgsea, GL_singscore, OM_singscore))

paired_test <- function(metric) {
  pr <- w[[paste0(metric, "_PR")]]; ln <- w[[paste0(metric, "_LN")]]
  delta <- ln - pr
  wt <- wilcox.test(ln, pr, paired = TRUE, exact = FALSE)
  # 配对 Cohen's d (差值均值/差值SD)
  d_eff <- mean(delta) / sd(delta)
  data.frame(metric = metric, mean_PR = mean(pr), mean_LN = mean(ln),
             mean_delta = mean(delta), median_delta = median(delta),
             n_up = sum(delta > 0), n_down = sum(delta < 0),
             cohen_d = d_eff, p_value = wt$p.value)
}
cat("\n=== 综合评分配对检验(n=8) ===\n")
res <- do.call(rbind, lapply(c("GL","OM"), paired_test))
print(res)

cat("\n=== 各算法子评分配对检验(稳健性) ===\n")
sub_metrics <- c("GL_ams","OM_ams","GL_aucell","OM_aucell",
                 "GL_ssgsea","OM_ssgsea","GL_singscore","OM_singscore")
res_sub <- do.call(rbind, lapply(sub_metrics, paired_test))
print(res_sub)
write.csv(res, file.path(ws, "results", "04_paired_tests_main.csv"), row.names = FALSE)
write.csv(res_sub, file.path(ws, "results", "04_paired_tests_by_algorithm.csv"),
          row.names = FALSE)

# ---- 象限在患者内 原发 vs 淋巴结 比例变化 ---------------------------------
cat("\n=== 患者级象限比例(原发 vs 淋巴结) ===\n")
qtab <- dp %>% group_by(lymph_node, quadrant) %>% summarise(cells = n(), .groups="drop") %>%
  group_by(lymph_node) %>% mutate(prop = cells / sum(cells))
print(as.data.frame(qtab))
write.csv(qtab, file.path(ws, "results", "04_quadrant_site_props.csv"), row.names = FALSE)

# 每患者 GL 高(GL_high)比例的配对
qgl <- dp %>% group_by(patient, lymph_node) %>%
  summarise(prop_GLhigh = mean(GL_score >= median(d$GL_score)), .groups="drop") %>%
  tidyr::pivot_wider(names_from = lymph_node, values_from = prop_GLhigh,
                     names_prefix = "site")
cat("\n各患者高GL比例: 原发 vs 淋巴结\n"); print(qgl)
p_gl <- wilcox.test(qgl$site1, qgl$site0, paired = TRUE, exact = FALSE)
cat("高GL比例配对Wilcoxon p =", p_gl$p.value,
    "| 平均增量:", mean(qgl$site1 - qgl$site0), "\n")
