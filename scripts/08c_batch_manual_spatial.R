# ============================================================================
# 08c_batch_manual_spatial.R
# 批量手动绘制 8 个 Visium 切片的 H&E + 代谢评分空间图
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

suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(tidyverse)})
ws <- root
sp_base <- file.path(hnscc_root, "data", "spatial_samples")
figdir <- file.path(ws, "figures")
samples <- list.files(sp_base)

plot_one <- function(sname) {
  s <- Load10X_Spatial(file.path(sp_base, sname),
                       filename = "filtered_feature_bc_matrix.h5",
                       image.name = "tissue_hires_image.png")
  imobj <- s@images[[1]]
  hs <- imobj@scale.factors$hires
  cent <- imobj@boundaries$centroids@coords
  spot <- data.frame(x = cent[,1] * hs, y = cent[,2] * hs)
  meta <- read.csv(file.path(ws, "results", "Spatial", paste0(sname, "_spots.csv")))
  stopifnot(nrow(spot) == nrow(meta))
  spot$HubScore <- meta$HubScore; spot$GL <- meta$GLspatial; spot$OM <- meta$OMspatial
  img <- imobj@image; h <- dim(img)[1]; w <- dim(img)[2]

  make <- function(value, title, colors) {
    ggplot(spot) +
      annotation_raster(img, 0, w, 0, h) +
      geom_point(aes(x, y, color = .data[[value]]), size = 1.3) +
      scale_color_gradientn(colors = colors) +
      scale_y_reverse() + coord_fixed() +
      labs(title = title, color = "") +
      theme_void(base_size = 11) +
      theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
  }
  p_hub <- make("HubScore", "Metabolic HubScore", c("#123B7A","white","#8A0612"))
  p_gl  <- make("GL", "Glycolysis/Lactate", c("white","#8A0612"))
  p_om  <- make("OM", "OXPHOS", c("white","#123B7A"))
  combined <- p_hub | p_gl | p_om
  ggsave(file.path(figdir, paste0("08_manual_", sname, ".pdf")),
         combined, width = 15, height = 6)
  ggsave(file.path(figdir, paste0("08_manual_", sname, ".png")),
         combined, width = 15, height = 6, dpi = 200)
  rm(s); gc()
}
for (sname in samples) { cat("plotting", sname, "\n"); plot_one(sname) }
cat("ALL DONE\n")
