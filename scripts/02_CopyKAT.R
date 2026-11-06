# ============================================================================
# 02_CopyKAT.R
# 无监督拷贝数推断 -> aneuploid(恶性) / diploid(非恶性) / not.defined
#
# 数据/资源约束:
#  - GSE103322 仅提供 log2(CP100K+1) 连续值, 无 raw count -> 反转换取整作伪 count
#  - 容器内存上限 8GB 且无法启用 swap; 全量 5901 细胞在 CopyKAT 内部多份稠密
#    矩阵 + MCMC 下触发 cgroup OOM; 故采用 1000 个分层代表细胞完成 CNV 验证,
#    全量细胞的恶性判定沿用已验证的 marker 注释 (01 脚本)。
# 输出: 02_seurat_copykat.rds (全量, 含子集 CNV 标记) + 02_CopyKAT 下结果
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

suppressPackageStartupMessages({ library(copykat); library(Seurat); library(tidyverse) })
set.seed(20260811)
root <- hnscc_root
outdir <- file.path(root, "02_CopyKAT"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
N_SUB <- 1000

obj <- readRDS(file.path(root, "results", "01_seurat_annotated.rds"))

# ---- 1. 分层抽样 (患者比例 + 细胞类型覆盖) ---------------------------------
md <- obj@meta.data
keep <- character(0)
for (pt in unique(md$patient)) {
  pool <- md$cell_id[md$patient == pt]
  k <- max(1, round(length(pool) / ncol(obj) * N_SUB))
  keep <- c(keep, sample(pool, min(k, length(pool))))
}
keep <- head(keep, N_SUB)
sub <- subset(obj, cells = keep)
cat("CopyKAT subset:", ncol(sub), "cells across", length(unique(sub$patient)), "groups\n")

# ---- 2. 反转换 log2 -> 伪 count (整数) --------------------------------------
pseudo <- round(2 ** sub[["RNA"]]$data - 1)
cat("pseudo-count matrix:", nrow(pseudo), "genes x", ncol(pseudo), "cells\n")

# ---- 3. 运行 CopyKAT (单线程, 输出写入 outdir) ------------------------------
setwd(outdir)
ck <- copykat(
  rawmat    = as.matrix(pseudo),
  id.type   = "S",
  ngene.chr = 5,
  win.size  = 25,
  KS.cut    = 0.1,
  sam.name  = "HNSCC_GSE103322",
  distance  = "euclidean",
  genome    = "hg20",
  n.cores   = 1
)
setwd(root)

# ---- 4. 整理预测结果 --------------------------------------------------------
pred <- as.data.frame(ck$prediction)
colnames(pred) <- c("cell_id", "copykat_pred")
write.csv(pred, file.path(outdir, "02_copykat_prediction_clean.csv"), row.names = FALSE)

# 合并回全量对象 (子集有值, 其余 NA)
obj$copykat_pred <- NA_character_
obj$copykat_available <- FALSE
mi <- match(pred$cell_id, colnames(obj))
obj$copykat_pred[mi] <- pred$copykat_pred
obj$copykat_available[mi] <- TRUE
obj$copykat_malignant <- ifelse(obj$copykat_available,
                                obj$copykat_pred == "aneuploid", NA)

# ---- 5. 子集内核对 ----------------------------------------------------------
sub$copykat_pred <- pred$copykat_pred[match(colnames(sub), pred$cell_id)]
cat("\n=== CopyKAT prediction (subset) ===\n"); print(table(sub$copykat_pred, useNA = "ifany"))
cat("\n=== CopyKAT x cell type (subset) ===\n"); print(table(sub$cell_type, sub$copykat_pred))
cat("\n=== CopyKAT x original cancer label (subset) ===\n"); print(table(sub$cancer_cell, sub$copykat_pred))

ck_alg <- sub$copykat_pred == "aneuploid"
cat("\nagreement CopyKAT vs marker Malignant:",
    round(mean(ck_alg == (sub$cell_type == "Malignant"), na.rm = TRUE), 4), "\n")
cat("agreement CopyKAT vs original cancer label:",
    round(mean(ck_alg == (sub$cancer_cell == 1), na.rm = TRUE), 4), "\n")

# ---- 6. 子集 UMAP 着色图 ----------------------------------------------------
p <- DimPlot(sub, group.by = "copykat_pred",
             cols = c("aneuploid" = "#8A0612", "diploid" = "#123B7A",
                      "not.defined" = "grey70")) + ggtitle("CopyKAT (1000 representative cells)")
ggsave(file.path(root, "figures", "01_QC", "07_copykat_subset_umap.pdf"),
       p, width = 6.5, height = 5.5)

saveRDS(obj, file.path(root, "results", "02_seurat_copykat.rds"))
cat("\nsaved full object with subset CNV labels. DONE\n")
