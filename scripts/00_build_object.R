# ============================================================================
# 00_build_object.R
# GSE103322 (Puram et al., Cell 2017) 口腔 HNSCC 单细胞 -> Seurat 对象
# 输入: data/GSE103322_HNSCC_all_data.txt.gz
# 数据性质: Smart-seq2/Fluidigm, 已归一化的 log2(CP100K+1) 连续表达值, 无 raw counts
# 输出: results/00_GSE103322_seurat_raw.rds
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
  library(Seurat); library(Matrix); library(tidyverse)
})
set.seed(20260811)

root <- hnscc_root
infile <- file.path(root, "data", "GSE103322_HNSCC_all_data.txt.gz")
outrds <- file.path(root, "results", "00_GSE103322_seurat_raw.rds")

# ---- 1. 读取完整矩阵 -------------------------------------------------------
raw <- read.delim(gzfile(infile), sep = "\t", check.names = FALSE,
                  stringsAsFactors = FALSE, row.names = 1)
message("raw dimension (rows x cells): ", nrow(raw), " x ", ncol(raw))

meta_labels <- c("processed by Maxima enzyme", "Lymph node",
                 "classified  as cancer cell", "classified as non-cancer cells",
                 "non-cancer cell type")
meta_df <- raw[meta_labels, , drop = FALSE]
expr_df <- raw[!rownames(raw) %in% meta_labels, , drop = FALSE]

# ---- 2. 清理基因名 (去单引号; 下划线统一为横线以匹配 Seurat 命名) ----------
genes <- str_remove_all(rownames(expr_df), "'")
genes <- gsub("_", "-", genes)
stopifnot(!any(duplicated(genes)))
rownames(expr_df) <- genes
X <- as.matrix(expr_df)
storage.mode(X) <- "double"
X <- Matrix(X, sparse = TRUE)
message("expression matrix (genes x cells): ", nrow(X), " x ", ncol(X))

# ---- 3. 构建规范细胞级元数据 -----------------------------------------------
cells <- colnames(X)
maxima <- as.character(meta_df["processed by Maxima enzyme", ])
lymph  <- as.integer(meta_df["Lymph node", ])
cancer <- as.integer(meta_df["classified  as cancer cell", ])
nctype <- as.character(meta_df["non-cancer cell type", ])

# 患者 ID: 规范化 HNSCC_17 -> HNSCC17; combo 单独标记
norm_patient <- function(cn) {
  if (startsWith(cn, "HNSCC_combo1")) return("combo1")
  if (startsWith(cn, "HNSCC_combo2")) return("combo2")
  if (startsWith(cn, "HNSCC_17")) return("HNSCC17")
  m <- regmatches(cn, regexpr("^[A-Z]+\\d+", cn))
  if (length(m) == 1) m else NA_character_
}
patient <- vapply(cells, norm_patient, character(1))
is_combo <- grepl("^combo", patient)

# 非癌细胞类型: 修正录入问题 (-Fibroblast -> Fibroblast); 0 表示癌细胞
celltype_raw <- nctype
celltype_raw[celltype_raw == "-Fibroblast"] <- "Fibroblast"
celltype_raw[celltype_raw == "0"] <- NA_character_

metadata <- data.frame(
  cell_id       = cells,
  patient       = patient,
  is_combo      = is_combo,
  cancer_cell   = cancer,                       # 1=原始注释癌细胞, 0=非癌
  lymph_node    = lymph,                        # 1=淋巴结, 0=原发
  maxima_enzyme = maxima,
  noncancer_type = celltype_raw,
  row.names     = cells,
  stringsAsFactors = FALSE
)
print(table(metadata$patient, useNA = "ifany"))
print(table(cancer, useNA = "ifany"))

# ---- 4. 创建 Seurat 对象 ---------------------------------------------------
# 数据已是 log2 归一化值: 同时作为 counts(占位) 与 data(实际 log 表达), 不再 NormalizeData
obj <- CreateSeuratObject(counts = X, meta.data = metadata, project = "HNSCC_GSE103322")
# 显式将 data layer 设为相同的 log2 矩阵 (避免对 log 值再次 log1p)
obj[["RNA"]]$data <- X
obj$log_transform <- "log2(CP100K+1), pre-normalized; no NormalizeData needed"
obj@misc$data_source <- "GSE103322 Puram Cell 2017; Smart-seq2; no raw counts"

# 基础 QC 指标 (基于 log 矩阵; 无 raw count 故 nFeature/nCount 来自现有矩阵)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^(RPS|RPL)")

message("final Seurat object: ", ncol(obj), " cells, ", nrow(obj), " genes")
saveRDS(obj, outrds)
message("saved: ", outrds)
