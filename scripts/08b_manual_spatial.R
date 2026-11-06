# ============================================================================
# 08b_manual_spatial.R
# 手动 ggplot 叠加 H&E 背景 + spot 代谢评分, 规避 Seurat VisiumV2 image 渲染问题
# 先在一个样本验证
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
base <- file.path(hnscc_root, "data", "spatial_samples", "P210325T1")
s <- Load10X_Spatial(base, filename = "filtered_feature_bc_matrix.h5",
                     image.name = "tissue_hires_image.png")
imobj <- s@images[[1]]
hires_scale <- imobj@scale.factors$hires
# full-res centroids
cent <- imobj@boundaries$centroids@coords
colnames(cent) <- c("x_full","y_full")
# 映射到 hires 坐标
spot <- data.frame(
  x = cent[, "x_full"] * hires_scale,
  y = cent[, "y_full"] * hires_scale)
# 读取已保存的该样本 spot 表(含 HubScore/GL/OM)
meta <- read.csv(file.path(res_dir, "Spatial", "P210325T1_spots.csv"))
stopifnot(nrow(spot) == nrow(meta))
spot$HubScore <- meta$HubScore; spot$GL <- meta$GLspatial; spot$OM <- meta$OMspatial

img <- imobj@image
h <- dim(img)[1]; w <- dim(img)[2]
cat("hires:", w, "x", h, "| spot x range:", range(spot$x), "y:", range(spot$y), "\n")

make <- function(value, title, colors) {
  ggplot(spot) +
    annotation_raster(img, xmin = 0, xmax = w, ymin = 0, ymax = h) +
    geom_point(aes(x = x, y = y, color = .data[[value]]), size = 1.4) +
    scale_color_gradientn(colors = colors) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = title, color = "") +
    theme_void(base_size = 11) +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
}
p_hub <- make("HubScore", "Metabolic HubScore", c("#123B7A","white","#8A0612"))
p_gl  <- make("GL", "GL spatial", c("white","#8A0612"))
p_om  <- make("OM", "OM spatial", c("white","#123B7A"))
combined <- (p_hub | p_gl | p_om)
ggsave("/tmp/manual_spatial_test.png", combined, width = 15, height = 6, dpi = 120)
cat("saved\n")
