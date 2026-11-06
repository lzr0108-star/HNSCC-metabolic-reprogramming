# ============================================================================
# 09_TCGA_metabolic.R
# TCGA-HNSC bulk 验证: 代谢共识基因 signature 预后 + 肿瘤/正常表达
# 内存友好: 逐个 fread STAR counts, 只取 hub 基因
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
  library(TCGAbiolinks); library(tidyverse); library(data.table)
  library(survival); library(survminer)
})
set.seed(20260812)
ws <- root
ann_root <- hnscc_root
base <- file.path(ann_root, "scripts", "GDCdata", "TCGA-HNSC",
                  "Transcriptome_Profiling", "Gene_Expression_Quantification")
dir.create(file.path(ws, "results", "TCGA"), showWarnings = FALSE, recursive = TRUE)
figdir <- file.path(ws, "figures")

q <- GDCquery(project = "TCGA-HNSC",
              data.category = "Transcriptome Profiling",
              data.type = "Gene Expression Quantification",
              workflow.type = "STAR - Counts")
map <- q$results[[1]] %>%
  select(id, case_id = cases.submitter_id, barcode = sample.submitter_id,
         sample_type) %>% distinct()
cat("files in query:", nrow(map), "\n")

cons <- read.csv(file.path(ws, "results", "ML", "consensus_genes.csv"))
hub <- head(cons$gene, 30)
ww_all <- setNames(cons$consensus_score, cons$gene)

read_one <- function(id) {
  f <- list.files(file.path(base, id), pattern = "\\.tsv$", full.names = TRUE)
  dt <- fread(f, select = c("gene_name","unstranded"), showProgress = FALSE)
  dt <- dt[gene_name != ""]
  lib <- sum(dt$unstranded)
  sub <- dt[gene_name %in% hub][, .(unstranded = sum(unstranded)), by = gene_name]
  list(lib = lib, sub = sub)
}
cnt_mat <- matrix(0, nrow = length(hub), ncol = nrow(map),
                  dimnames = list(hub, map$barcode))
libsize <- setNames(numeric(nrow(map)), map$barcode)
for (i in seq_len(nrow(map))) {
  r <- read_one(map$id[i]); libsize[i] <- r$lib
  cnt_mat[r$sub$gene_name, map$barcode[i]] <- r$sub$unstranded
  if (i %% 100 == 0) cat("read", i, "of", nrow(map), "\n")
}
cat("count matrix:", nrow(cnt_mat), "x", ncol(cnt_mat), "\n")
logcpm <- log2(sweep(cnt_mat, 2, libsize, "/") * 1e6 + 1)
hub <- intersect(hub, rownames(logcpm))
# 剔除无变异(标准差为0, 多为全0未检出)基因, 避免 scale 产生 NA 污染 signature
gene_sd <- apply(logcpm[hub, , drop = FALSE], 1, sd)
drop_g <- hub[gene_sd == 0 | is.na(gene_sd)]
if (length(drop_g) > 0) cat("剔除无变异基因:", paste(drop_g, collapse=", "), "\n")
hub <- hub[!(hub %in% drop_g)]; ww <- ww_all[hub]
cat("用于 signature 的基因数:", length(hub), "\n")

clinic <- tryCatch(GDCquery_clinic(project = "TCGA-HNSC", type = "clinical"),
                   error = function(e) NULL)
tries <- 0
while ((is.null(clinic) || nrow(clinic) < 100) && tries < 5) {
  tries <- tries + 1; cat("clinic 获取异常, 重试", tries, "\n"); Sys.sleep(5)
  clinic <- tryCatch(GDCquery_clinic(project = "TCGA-HNSC", type = "clinical"),
                     error = function(e) NULL)
}
stopifnot(!is.null(clinic), nrow(clinic) >= 100)
cd <- map %>% left_join(clinic, by = c("case_id" = "submitter_id")) %>%
  mutate(time = as.numeric(ifelse(!is.na(days_to_death), days_to_death,
                                  days_to_last_follow_up)),
         status = ifelse(vital_status == "Dead", 1, 0),
         stage = ajcc_pathologic_stage, grade = tumor_grade,
         age = age_at_index, gender = sex_at_birth,
         is_tumor = sample_type == "Primary Tumor")
rownames(cd) <- cd$barcode

Z <- scale(t(logcpm[hub, , drop = FALSE]))
cat("DIAG Z dim:", dim(Z), "| hub len:", length(hub),
    "| 列名匹配:", all(colnames(Z) == hub), "\n")
cat("DIAG Z NA 数:", sum(is.na(Z)), "\n")
sig <- as.numeric(Z[, hub, drop = FALSE] %*% ww[hub] / sum(ww[hub]))
cat("DIAG sig length:", length(sig), "NA:", sum(is.na(sig)),
    "| cd nrow:", nrow(cd), "\n")
cd$signature <- as.numeric(scale(sig))
cd$sig_group <- ifelse(cd$signature >= median(cd$signature, na.rm = TRUE), "High","Low")
cat("DIAG signature NA:", sum(is.na(cd$signature)), "\n")

cat("DIAG nrow(cd):", nrow(cd), "| is_tumor:", sum(cd$is_tumor),
    "| time非空:", sum(!is.na(cd$time)), "| status非空:", sum(!is.na(cd$status)), "\n")
cat("DIAG vital_status:\n"); print(table(cd$vital_status, useNA = "ifany"))
surv_df <- cd %>% filter(is_tumor, !is.na(time), !is.na(status))
cat("DIAG surv_df nrow:", nrow(surv_df), "\n")
cat("DIAG time class:", class(surv_df$time), "range:", range(surv_df$time), "\n")
cat("DIAG status class:", class(surv_df$status), "\n"); print(table(surv_df$status, useNA="ifany"))
cat("DIAG sig_group:\n"); print(table(surv_df$sig_group, useNA="ifany"))
cat("DIAG Surv nonmissing:", sum(!is.na(surv_df$time) & !is.na(surv_df$status)), "\n")
fit_km <- survfit(Surv(time, status) ~ sig_group, data = surv_df)
lr <- surv_pvalue(fit_km, data = surv_df)
cat("KM log-rank p:", lr$pval, "\n")
p1 <- ggsurvplot(fit_km, data = surv_df, pval = TRUE, risk.table = TRUE,
                 palette = c("#8A0612","#123B7A"),
                 title = "TCGA-HNSC metabolic signature KM")
p_km <- cowplot::plot_grid(p1$plot, p1$table, ncol = 1,
                           rel_heights = c(3,1), align = "v")
ggsave(file.path(figdir, "09_KM_signature.pdf"), p_km, width = 8, height = 9)
ggsave(file.path(figdir, "09_KM_signature.png"), p_km, width = 8, height = 9, dpi = 200)

cox1 <- coxph(Surv(time, status) ~ signature + age + gender, data = surv_df)
cox_tab <- as.data.frame(broom::tidy(cox1, exponentiate = TRUE, conf.int = TRUE))
write.csv(cox_tab, file.path(ws, "results", "TCGA", "cox_regression.csv"), row.names = FALSE)
print(cox_tab)

single <- lapply(head(hub, 12), function(g) {
  d <- surv_df; d$g_expr <- as.numeric(logcpm[g, d$barcode])
  d$g_group <- ifelse(d$g_expr >= median(d$g_expr), "High","Low")
  f <- survfit(Surv(time, status) ~ g_group, data = d)
  data.frame(gene = g, logrank_p = surv_pvalue(f, data = d)$pval)
})
single <- bind_rows(single); single$fdr <- p.adjust(single$logrank_p, "BH")
write.csv(single, file.path(ws, "results", "TCGA", "single_gene_survival.csv"), row.names = FALSE)

tn <- cd %>% filter(sample_type %in% c("Primary Tumor","Solid Tissue Normal"))
tn_res <- lapply(hub, function(g) {
  vals <- data.frame(g = as.numeric(logcpm[g, tn$barcode]),
                     type = ifelse(tn$sample_type == "Primary Tumor","Tumor","Normal"))
  data.frame(gene = g, tumor_mean = mean(vals$g[vals$type=="Tumor"]),
             normal_mean = mean(vals$g[vals$type=="Normal"]),
             log2FC = mean(vals$g[vals$type=="Tumor"]) - mean(vals$g[vals$type=="Normal"]),
             pval = suppressWarnings(wilcox.test(g ~ type, vals)$p.value))
})
tn_res <- bind_rows(tn_res); tn_res$fdr <- p.adjust(tn_res$pval, "BH")
write.csv(tn_res, file.path(ws, "results", "TCGA", "tumor_normal_DE.csv"), row.names = FALSE)
cat("\n=== tumor vs normal (前10) ===\n"); print(head(tn_res,10))

p2 <- tn_res %>% ggplot(aes(reorder(gene, log2FC), log2FC, fill = log2FC > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("#123B7A","#8A0612")) +
  labs(x = NULL, y = "Tumor vs Normal log2 difference") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(figdir, "09_tumor_normal_bar.pdf"), p2, width = 8, height = 8)
ggsave(file.path(figdir, "09_tumor_normal_bar.png"), p2, width = 8, height = 8, dpi = 200)

p3 <- surv_df %>% filter(!is.na(stage)) %>%
  ggplot(aes(stage, signature)) +
  geom_boxplot(fill = "#123B7A", color = "#123B7A", alpha = 0.85, outlier.size = 0.6) +
  coord_flip() + labs(x = NULL, y = "Metabolic signature") + theme_bw()
ggsave(file.path(figdir, "09_signature_by_stage.pdf"), p3, width = 8, height = 6)

write.csv(surv_df %>% select(barcode, case_id, time, status, stage, grade,
                             signature, sig_group),
          file.path(ws, "results", "TCGA", "patient_signature.csv"), row.names = FALSE)
saveRDS(list(km_p = lr$pval, cox = cox_tab, single = single, tn = tn_res),
        file.path(ws, "results", "TCGA", "TCGA_summary.rds"))
cat("\nDONE\n")
