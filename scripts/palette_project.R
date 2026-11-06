# ============================================================================
# palette_project.R
# 项目专用调色板: 自定义、色盲友好, 区别于 ggsci NPG 与 Seurat/ggplot 默认配色
# 由 01_QC_cluster.R / 01b_annotate.R / 02_CopyKAT.R 及重新着色脚本统一 source
# ============================================================================

suppressPackageStartupMessages(library(colorspace))
suppressPackageStartupMessages(library(ggplot2))

# 11 大细胞类型 (命名定性色)
cell_pal <- c(
  "Malignant"   = "#A8325A",  # 深莓红
  "CAF"         = "#E0A852",  # 焦糖金
  "T/NK"        = "#2F6F8E",  # 深青蓝
  "Treg"        = "#6A51A3",  # 深紫
  "B cell"      = "#1F9E89",  # 翡翠松石绿
  "Plasma"      = "#394F8A",  # 靛蓝
  "Endothelial" = "#CC5A33",  # 赤陶橙
  "Macrophage"  = "#C06C84",  # 玫瑰粉
  "Mast"        = "#8C6D46",  # 棕褐
  "Dendritic"   = "#5B8AB0",  # 雾蓝
  "Myocyte"     = "#7CA982"   # 橄榄绿
)

# 二元原始 cancer 标签 (0 非恶性 / 1 恶性)
cancer_pal <- c("0" = "#4A90A4", "1" = "#C1476B")

# CopyKAT 三类
copykat_pal <- c("aneuploid" = "#A8325A", "diploid" = "#2F6F8E",
                 "not.defined" = "#B5B5B5")

# 多类别 (patient / cluster / subtype): 固定 HCL 参数, 可复现
multi_pal <- function(n) {
  colorspace::qualitative_hcl(n, h = c(0, 300), c = 62, l = 60)
}

# marker 点图连续色阶 (浅灰 -> 深青)
dot_scale <- function(...) {
  ggplot2::scale_color_gradient(low = "#E6E6E6", high = "#2F6F8E", ...)
}
