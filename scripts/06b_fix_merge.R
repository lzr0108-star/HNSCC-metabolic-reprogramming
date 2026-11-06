# ============================================================================
# 06b_fix_merge.R
# 重新加载并合并 HighGL/LowGL CellChat, 保存 merge, 稳健导出通讯表
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

suppressPackageStartupMessages({library(CellChat); library(tidyverse)})
ws <- root
cc_high <- readRDS(file.path(ws, "results", "CellChat", "cellchat_HighGL.rds"))
cc_low  <- readRDS(file.path(ws, "results", "CellChat", "cellchat_LowGL.rds"))
cc_list <- list(HighGL = cc_high, LowGL = cc_low)
cc_merge <- mergeCellChat(cc_list, add.names = c("HighGL","LowGL"), cell.prefix = TRUE)
saveRDS(cc_merge, file.path(ws, "results", "CellChat", "cellchat_merged.rds"))
cat("merge对象已保存\n")

# 分别导出通讯表后合并
ch <- subsetCommunication(cc_high) %>% mutate(dataset = "HighGL")
cl <- subsetCommunication(cc_low)  %>% mutate(dataset = "LowGL")
cat("HighGL通讯条目:", nrow(ch), "| LowGL:", nrow(cl), "\n")
common_cols <- intersect(colnames(ch), colnames(cl))
merged_comm <- bind_rows(ch[, common_cols], cl[, common_cols])
write.csv(merged_comm, file.path(ws, "results", "CellChat", "communication_merged.csv"),
          row.names = FALSE)
cat("合并通讯表已保存:", nrow(merged_comm), "行\n")

# 通路级总强度对比
path_h <- cc_high@netP$pathways; path_l <- cc_low@netP$pathways
cat("\nHighGL通路数:", length(path_h), "| LowGL:", length(path_l), "\n")
# 每个通路总强度(netP prob 求和)
path_strength <- function(cc, name) {
  pw <- cc@netP$prob
  data.frame(pathway = dimnames(pw)[[3]],
             strength = apply(pw, 3, sum), dataset = name)
}
ps <- bind_rows(path_strength(cc_high, "HighGL"), path_strength(cc_low, "LowGL"))
w <- ps %>% pivot_wider(names_from = dataset, values_from = strength,
                        values_fill = 0) %>%
  mutate(delta = HighGL - LowGL) %>% arrange(desc(abs(delta)))
write.csv(w, file.path(ws, "results", "CellChat", "pathway_strength_compare.csv"),
          row.names = FALSE)
cat("\n=== HighGL 相对 LowGL 增强 Top10 通路 ===\n")
print(head(w[order(-w$delta), ], 10))
cat("\n=== HighGL 相对 LowGL 减弱 Top10 通路 ===\n")
print(head(w[order(w$delta), ], 10))
