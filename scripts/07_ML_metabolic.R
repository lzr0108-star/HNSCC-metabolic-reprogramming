# ============================================================================
# 07_ML_metabolic.R
# 10 机器学习 + 6 特征选择 -> 高/低糖酵解共识预测基因
# 关键红线: 特征池剔除全部 GL/OM 评分基因, 防止循环论证
# 按 patient 7:3 拆分; 二分类 HighGL vs LowGL
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
  library(Seurat); library(caret); library(tidyverse); library(pROC); library(glmnet)
})
set.seed(20260811)
ws <- root
dir.create(file.path(ws, "results", "ML"), showWarnings = FALSE, recursive = TRUE)
figdir <- file.path(ws, "figures")

z <- readRDS(file.path(ws, "results", "03_malignant_grouped.rds"))
mal <- z$mal
gs <- readRDS(file.path(ws, "results", "03_metabolic_genesets.rds"))
exclude_genes <- unique(c(gs$GL, gs$OM))
cat("需剔除评分基因数:", length(exclude_genes), "\n")

combo_col <- if ("is_combo" %in% colnames(mal@meta.data)) "is_combo" else NULL
keep <- mal$GL_group %in% c("HighGL","LowGL")
if (!is.null(combo_col)) keep <- keep & !mal$is_combo
seu <- subset(mal, cells = colnames(mal)[keep])
y <- factor(ifelse(seu$GL_group == "HighGL", "High", "Low"), levels = c("Low","High"))
names(y) <- colnames(seu)
cat("analysis cells:", length(y), "| High:", sum(y=="High"), "Low:", sum(y=="Low"), "\n")

X <- t(as.matrix(seu[["RNA"]]$data))
X <- X[, colSums(X) > 0, drop = FALSE]
# === 红线: 剔除评分基因 ===
before <- ncol(X)
X <- X[, !(colnames(X) %in% exclude_genes), drop = FALSE]
cat("剔除评分基因:", before - ncol(X), "| 剩余候选:", ncol(X), "\n")
clean <- make.names(colnames(X))
name_map <- setNames(colnames(X), clean)
colnames(X) <- clean
pat <- seu$patient
rm(seu, mal); gc()

# ---- patient 7:3 拆分 -------------------------------------------------------
dom <- tapply(as.numeric(y) - 1, pat, mean)
high_pats <- names(dom)[dom >= 0.5]; low_pats <- names(dom)[dom < 0.5]
split_by_count <- function(pgroup, ratio = 0.7) {
  pc <- sort(table(pat[pat %in% pgroup]), decreasing = TRUE)
  target <- sum(pc) * ratio; tr <- c(); acc <- 0
  for (pn in names(pc)) if (acc < target) { tr <- c(tr, pn); acc <- acc + pc[[pn]] }
  tr
}
train_pats <- c(split_by_count(high_pats), split_by_count(low_pats))
test_pats <- setdiff(c(high_pats, low_pats), train_pats)
tr_idx <- names(pat)[pat %in% train_pats]; te_idx <- names(pat)[pat %in% test_pats]
cat("train patients:", length(train_pats), "| test:", length(test_pats), "\n")
print(table(y[tr_idx])); print(table(y[te_idx]))

# ---- 6 特征选择 (train only) ------------------------------------------------
topn <- 60
Xtr <- X[tr_idx, , drop = FALSE]; ytr <- y[tr_idx]; ynum <- as.numeric(ytr) - 1
mk <- sapply(colnames(Xtr), function(g)
  suppressWarnings(wilcox.test(Xtr[ytr=="High",g], Xtr[ytr=="Low",g])$p.value))
lfc <- colMeans(Xtr[ytr=="High",,drop=FALSE]) - colMeans(Xtr[ytr=="Low",,drop=FALSE])
fs_wilcox <- names(sort(mk[p.adjust(mk,"BH") < 0.05 & abs(lfc) >= 0.5]))[1:topn]
set.seed(20260811); cv_fold <- sample(rep(1:5, length.out=nrow(Xtr)))
cv_l <- cv.glmnet(Xtr, ynum, family="binomial", alpha=1, foldid=cv_fold)
cf <- as.matrix(coef(cv_l, s="lambda.min")); fs_lasso <- rownames(cf)[cf[,1]!=0][-1]
cv_e <- cv.glmnet(Xtr, ynum, family="binomial", alpha=0.5, foldid=cv_fold)
cfe <- as.matrix(coef(cv_e, s="lambda.min")); fs_enet <- rownames(cfe)[cfe[,1]!=0][-1]
rfm <- randomForest::randomForest(Xtr, ytr, importance=TRUE, ntree=500)
imp <- randomForest::importance(rfm, type=1)[,1]; fs_rf <- names(sort(imp, decreasing=TRUE))[1:topn]
cors <- apply(Xtr, 2, function(g) suppressWarnings(cor(g, ynum)))
fs_cor <- names(sort(abs(cors), decreasing=TRUE))[1:topn]
rfe_genes <- head(unique(c(fs_wilcox, fs_lasso, fs_rf, fs_cor)), 200)
rfe_ctrl <- rfeControl(functions=caretFuncs, method="cv", number=3, verbose=FALSE)
rfe_res <- rfe(Xtr[, rfe_genes, drop=FALSE], ytr, sizes=c(20,40,60,80),
               rfeControl=rfe_ctrl, method="svmRadial", fit=FALSE)
fs_svmrfe <- predictors(rfe_res)
feature_sets <- list(Wilcox=fs_wilcox, LASSO=fs_lasso, ElasticNet=fs_enet,
                     RF=fs_rf, Correlation=fs_cor, SVMRFE=fs_svmrfe)
feature_sets <- lapply(feature_sets, function(v) v[!is.na(v) & v!=""])
print(lengths(feature_sets))
pool <- unique(unlist(feature_sets)); cat("union pool:", length(pool), "\n")

# ---- 10 算法训练 ------------------------------------------------------------
methods10 <- c("rf","glmnet","svmRadial","pls","knn","rpart","nnet","gbm","LogitBoost","nb")
tr_ctrl <- trainControl(method="cv", number=5, classProbs=TRUE,
                        summaryFunction=twoClassSummary, savePredictions="final")
Xtr_p <- Xtr[, pool, drop=FALSE]; Xte_p <- X[te_idx, pool, drop=FALSE]
models <- list(); pred_prob <- list()
for (mth in methods10) {
  cat("training:", mth, "\n")
  fit <- tryCatch(suppressWarnings(train(Xtr_p, ytr, method=mth, metric="ROC",
                                         trControl=tr_ctrl, tuneLength=3)),
                  error=function(e){cat(" ERROR:", conditionMessage(e),"\n"); NULL})
  if (!is.null(fit)) { models[[mth]] <- fit
    pred_prob[[mth]] <- as.numeric(suppressWarnings(predict(fit, Xte_p, type="prob")[,"High"])) }
}
ytest <- as.numeric(y[te_idx]) - 1
eval_df <- do.call(rbind, lapply(names(pred_prob), function(mth){
  p <- pred_prob[[mth]]; pc <- ifelse(p>=0.5,1,0)
  data.frame(method=mth, AUC=as.numeric(auc(roc(ytest,p,quiet=TRUE))),
             Accuracy=mean(pc==ytest), Sensitivity=mean(pc[ytest==1]==1),
             Specificity=mean(pc[ytest==0]==0)) })) %>% arrange(-AUC)
write.csv(eval_df, file.path(ws,"results/ML/model_performance.csv"), row.names=FALSE)
print(eval_df, row.names=FALSE)

# ---- 共识基因 ---------------------------------------------------------------
fs_count <- sapply(pool, function(g) sum(sapply(feature_sets, function(s) g %in% s)))
var_imp <- sapply(names(models), function(mth){
  v <- tryCatch(varImp(models[[mth]])$importance, error=function(e) NULL)
  if (is.null(v) || !"High" %in% colnames(v)) return(setNames(rep(0,length(pool)),pool))
  z <- v[pool,"High"]; z[is.na(z)] <- 0
  if (max(z,na.rm=TRUE)>0) z/max(z,na.rm=TRUE) else z })
imp_mean <- rowMeans(as.matrix(var_imp), na.rm=TRUE)
consensus <- data.frame(gene=name_map[pool], feature_selection_count=fs_count,
  mean_variable_importance=round(imp_mean,4),
  consensus_score=round(fs_count/6*0.5 + imp_mean*0.5,4)) %>% arrange(-consensus_score)
write.csv(consensus, file.path(ws,"results/ML/consensus_genes.csv"), row.names=FALSE)
cat("Top 20 consensus (均为非评分基因):\n"); print(head(consensus,20), row.names=FALSE)

# ---- 出图 -------------------------------------------------------------------
p1 <- eval_df %>% pivot_longer(AUC:Specificity, names_to="metric", values_to="value") %>%
  ggplot(aes(reorder(method,value), value, fill=method)) + geom_col() + coord_flip() +
  facet_wrap(~metric) + theme_bw() + theme(legend.position="none")
ggsave(file.path(figdir,"07_performance.pdf"), p1, width=9, height=6)
ggsave(file.path(figdir,"07_performance.png"), p1, width=9, height=6, dpi=300)
p3 <- head(consensus,25) %>%
  ggplot(aes(reorder(gene,consensus_score), consensus_score, fill=feature_selection_count)) +
  geom_col() + coord_flip() + scale_fill_gradient(low="#D8B96A", high="#8A0612") +
  labs(x=NULL, y="Consensus score") + theme_bw()
ggsave(file.path(figdir,"07_consensus_score.pdf"), p3, width=7, height=7)
ggsave(file.path(figdir,"07_consensus_score.png"), p3, width=7, height=7, dpi=300)

saveRDS(list(performance=eval_df, consensus=consensus, feature_sets=feature_sets,
             train_patients=train_pats, test_patients=test_pats),
        file.path(ws,"results/ML/ML_summary.rds"))
cat("\nDONE\n")
