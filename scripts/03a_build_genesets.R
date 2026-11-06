# ============================================================================
# 03a_build_genesets.R
# 构建 GL(糖酵解/乳酸) 与 OM(线粒体 OXPHOS) 双轴基因集
# 来源: Hallmark / KEGG / Reactome (msigdbr) + 手工乳酸转运/乳酸化基因
# 输出: results/03_metabolic_genesets.rds
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

suppressPackageStartupMessages({library(Seurat); library(dplyr); library(msigdbr)})
set.seed(20260811)
ws <- root
obj <- readRDS(file.path(hnscc_root, "results", "01_seurat_annotated.rds"))
available_genes <- rownames(obj)
cat("对象基因数:", length(available_genes), "\n")

H  <- msigdbr(species = "Homo sapiens", collection = "H")
C2 <- msigdbr(species = "Homo sapiens", collection = "C2")
getset <- function(db, name) db %>% filter(gs_name == name) %>% pull(gene_symbol) %>% unique()

# ---- GL 轴各独立子集合 -----------------------------------------------------
gl_sub <- list(
  Hallmark_glycolysis      = getset(H,  "HALLMARK_GLYCOLYSIS"),
  KEGG_glycolysis          = getset(C2, "KEGG_GLYCOLYSIS_GLUCONEOGENESIS"),
  KEGG_pyruvate            = getset(C2, "KEGG_PYRUVATE_METABOLISM"),
  Reactome_glycolysis_reg  = getset(C2, "REACTOME_REGULATION_OF_GLYCOLYSIS_BY_FRUCTOSE_2_6_BISPHOSPHATE_METABOLISM"),
  Reactome_pyruvate        = getset(C2, "REACTOME_PYRUVATE_METABOLISM"),
  Lactate_manual           = c("LDHA","LDHB","SLC16A1","SLC16A3","SLC16A7","EP300")
)
# ---- OM 轴各独立子集合 -----------------------------------------------------
om_sub <- list(
  Hallmark_oxphos          = getset(H,  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
  KEGG_oxphos              = getset(C2, "KEGG_OXIDATIVE_PHOSPHORYLATION"),
  KEGG_tca                 = getset(C2, "KEGG_CITRATE_CYCLE_TCA_CYCLE"),
  Reactome_tca             = getset(C2, "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE"),
  Reactome_resp_chain      = getset(C2, "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"),
  Reactome_tca_maturation  = getset(C2, "REACTOME_MATURATION_OF_TCA_ENZYMES_AND_REGULATION_OF_TCA_CYCLE")
)

# ---- 过滤: 去 MT 编码基因, 仅保留对象中存在的基因 ---------------------------
clean <- function(g) {
  g <- unique(g)
  g <- g[!grepl("^MT-", g)]
  g[g %in% available_genes]
}
gl_sub <- lapply(gl_sub, clean)
om_sub <- lapply(om_sub, clean)

GL <- unique(unlist(gl_sub, use.names = FALSE))
OM <- unique(unlist(om_sub, use.names = FALSE))
cat("\n=== GL 子集合(过滤后基因数) ===\n"); print(sapply(gl_sub, length))
cat("GL 综合基因数:", length(GL), "\n")
cat("\n=== OM 子集合(过滤后基因数) ===\n"); print(sapply(om_sub, length))
cat("OM 综合基因数:", length(OM), "\n")
cat("\n两轴重叠基因数:", length(intersect(GL, OM)), "\n")
print(intersect(GL, OM))

genesets <- list(GL = GL, OM = OM, gl_sub = gl_sub, om_sub = om_sub)
saveRDS(genesets, file.path(ws, "results", "03_metabolic_genesets.rds"))
cat("\n已保存基因集\n")
