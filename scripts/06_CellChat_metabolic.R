# ============================================================================
# 06_CellChat_metabolic.R
# 高糖酵解(HighGL) vs 低糖酵解(LowGL) 恶性细胞通讯比较
# HighGL 数据集: HighGL 恶性 + 全部非恶性; LowGL 同理
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
  library(Seurat); library(CellChat); library(patchwork); library(tidyverse)
})
set.seed(20260811)
ws <- root
ann_root <- hnscc_root
dir.create(file.path(ws, "results", "CellChat"), showWarnings = FALSE, recursive = TRUE)
figdir <- file.path(ws, "figures"); dir.create(figdir, showWarnings = FALSE)

ld <- new.env(); lazyLoad(file.path(system.file(package = "CellChat"), "data", "Rdata"),
                         envir = ld)
CellChatDB <- ld$CellChatDB.human

# 预读分组并映射
z <- readRDS(file.path(ws, "results", "03_malignant_grouped.rds"))
gl_map <- setNames(as.character(z$mal$GL_group), colnames(z$mal))

prepare <- function() {
  src <- readRDS(file.path(ann_root, "results", "01_seurat_annotated.rds"))
  src[["RNA"]]$scale.data <- NULL
  src$GL_group <- NA
  src$GL_group[names(gl_map)] <- gl_map
  src$cc_label <- ifelse(src$cell_type == "Malignant", "Malignant",
                         as.character(src$cell_type))
  src
}
src <- prepare()
high_cells <- colnames(src)[src$cc_label != "Malignant" | src$GL_group == "HighGL"]
low_cells  <- colnames(src)[src$cc_label != "Malignant" | src$GL_group == "LowGL"]
cat("HighGL dataset:", length(high_cells), "| LowGL dataset:", length(low_cells), "\n")
rm(src); gc()

run_cellchat <- function(cells, label) {
  seu <- subset(prepare(), cells = cells); gc()
  cc <- createCellChat(object = seu, group.by = "cc_label", assay = "RNA")
  cc@DB <- CellChatDB
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1,
                          population.size = FALSE, raw.use = TRUE, seed.use = 20260811)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  saveRDS(cc, file.path(ws, "results", "CellChat", paste0("cellchat_", label, ".rds")))
  cc
}
cc_high <- run_cellchat(high_cells, "HighGL")
cc_low  <- run_cellchat(low_cells, "LowGL")

cc_list <- list(HighGL = cc_high, LowGL = cc_low)
cc_merge <- mergeCellChat(cc_list, add.names = c("HighGL","LowGL"), cell.prefix = TRUE)

# ---- 出图 ------------------------------------------------------------------
for (nm in c("HighGL","LowGL")) {
  cc <- cc_list[[nm]]
  pdf(file.path(figdir, paste0("06_circle_", nm, ".pdf")), width = 7, height = 7)
  try(netVisual_circle(cc@net$count, vertex.weight = rowSums(cc@net$count),
                       weight.scale = TRUE, title.name = paste(nm, "counts")))
  try(netVisual_circle(cc@net$weight, vertex.weight = rowSums(cc@net$weight),
                       weight.scale = TRUE, title.name = paste(nm, "weight")))
  invisible(dev.off())
}
try(ggsave(file.path(figdir, "06_interaction_count.pdf"),
           compareInteractions(cc_merge, show.legend = FALSE), width = 5, height = 4))
try(ggsave(file.path(figdir, "06_diff_interaction.pdf"),
           netVisual_diffInteraction(cc_merge, comparison = c(1,2), weight.scale = TRUE),
           width = 6, height = 5))
for (nm in c("HighGL","LowGL")) {
  pdf(file.path(figdir, paste0("06_role_heatmap_", nm, ".pdf")), width = 8, height = 6)
  try(netAnalysis_signalingRole_heatmap(cc_list[[nm]], pattern = "outgoing"))
  try(netAnalysis_signalingRole_heatmap(cc_list[[nm]], pattern = "incoming"))
  invisible(dev.off())
}
try(ggsave(file.path(figdir, "06_rankNet.pdf"),
           rankNet(cc_merge, mode = "comparison", comparison = c(1,2),
                   stacked = TRUE, do.stat = TRUE), width = 8, height = 8))
# bubble: Malignant 作为发送者(动态索引)
mal_idx <- which(levels(cc_high@idents) == "Malignant")
common_path <- intersect(cc_high@netP$pathways, cc_low@netP$pathways)
try(ggsave(file.path(figdir, "06_bubble_malignant.pdf"),
           netVisual_bubble(cc_merge, sources.use = mal_idx, targets.use = NULL,
                            comparison = c(1,2), signaling = head(common_path, 12),
                            remove.isolate = FALSE), width = 9, height = 7))

comm <- subsetCommunication(cc_merge)
write.csv(comm, file.path(ws, "results", "CellChat", "communication_merged.csv"),
          row.names = FALSE)
saveRDS(cc_merge, file.path(ws, "results", "CellChat", "cellchat_merged.rds"))
cat("\nDONE\n")
