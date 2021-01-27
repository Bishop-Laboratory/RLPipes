## Original model building script ##

# For reference, this is the script I used initially to preprocess the data,
# assign labels, and build the classifier and regression quality models. 
# There are several issues with my approach that are addressed in the modules.
# Firstly, I initially allowed only human data and only DRIP-family data. 
# However, I think that we should be able to include more datasets/species with some tuning. 
# Secondly, ML tasks like feature selection, feature engineering, model selection, model validation,
# and parameter tuning were not done in a proper manner -- so this also needs to be addressed. 

# Scores
library(tidyverse)
library(ChIPpeakAnno)
library(factoextra)
rmap_samples <- read_csv("misc/rmap_full_11_25.csv")
rmap_features <- rmap_samples %>%
  mutate(has_control = ifelse(is.na(control), 0, 1)) %>%
  dplyr::select(clean_name, # corr_median, corr_skewness, corr_kurtosis, 
                rlfs_pval,
                rlfs_median, rlfs_skewness, rlfs_kurtosis, `3UTR__LogP enrichment (+values depleted)`, `Retroposon__LogP enrichment (+values depleted)`,
                `miRNA__LogP enrichment (+values depleted)`, `TTS__LogP enrichment (+values depleted)`, `SINE__LogP enrichment (+values depleted)`,
                `LINE__LogP enrichment (+values depleted)`, `tRNA__LogP enrichment (+values depleted)`, `rRNA__LogP enrichment (+values depleted)`,
                `Exon__LogP enrichment (+values depleted)`, `Intron__LogP enrichment (+values depleted)`, `Intergenic__LogP enrichment (+values depleted)`,
                `Promoter__LogP enrichment (+values depleted)`, `5UTR__LogP enrichment (+values depleted)`,
                `Satellite__LogP enrichment (+values depleted)`, `MACS2__Percentage Peaks Overlapping`,
                `EPIC2__Percentage Peaks Overlapping`,
                gc_content, pct_aligned,
                percent_passing, pct_duplicated) %>%
  mutate_at(.vars = vars(contains('LogP')), .funs = function(x) {as.numeric(as.character(x)) * -1})
colnames(rmap_features) <- gsub(colnames(rmap_features), pattern = "__L.+", replacement = "")
colnames(rmap_features) <- gsub(colnames(rmap_features), pattern = "__P.+", replacement = " % Overlap")
colnames(rmap_features) <- gsub(colnames(rmap_features), pattern = "rlfs_", replacement = "RLFS ")
colnames(rmap_features) <- gsub(colnames(rmap_features), pattern = "corr_", replacement = "Correlation ")
colnames(rmap_features) <- gsub(colnames(rmap_features), pattern = "gc_", replacement = "GC ")
rmap_features <- rmap_featuresfull %>%
  column_to_rownames(var = "clean_name") %>%
  as.matrix()
rmap_features <- apply(rmap_features, 1:2, as.numeric)
rmap_features <- apply(rmap_features, 1:2, function(x) {ifelse(is.na(x), 0, x)})
rmap_features <- rmap_features[which(rowSums(rmap_features) != 0), which(colSums(rmap_features) != 0)]

feat_cor <- cor(rmap_features)
feat_cor <- apply(feat_cor, 1:2, function(x) {ifelse(is.na(x), 0, x)})
pheatmap::pheatmap(feat_cor, filename = "results/quality_feature_correlation.png",
                   height = 8, width = 10)
dev.off()
rmap_featuresscale <- scale(rmap_features)
rmap_featuresscale <- apply(rmap_featuresscale, 1:2, function(x) {ifelse(is.na(x), 0, x)})
rmap_featuresscale <- rmap_featuresscale[, which(colSums(rmap_featuresscale) != 0)]

pcafeat <- prcomp(rmap_featuresscale)

scree <- fviz_eig(pcafeat) +
  theme_bw(base_size = 15) +
  xlab("Principal Component")
scree + ggsave(filename = "results/pca_scree_plot.png")

var_exp <- scree$data
clusts <- hclust(dist(rmap_featuresscale), method = 'single')


png("results/dendrogram_rmap.png", height = 10, width = 50, units = "in", res = 300)
dend <- as.dendrogram(clusts)
par(mar = c(20, 5, 5,5))
plot(dend,  xlab = "Height",
     horiz = F)
dev.off()
dev.off()

clusters <- cutree(clusts,  k = 40)
pct_pc1 <- round(var_exp$eig[1], 2)
pct_pc2 <- round(var_exp$eig[2], 2)
rmap_clust <- rmap_featuresfull %>%
  mutate(PC1 = pcafeat$x[,c(1)]) %>%
  mutate(PC2 = pcafeat$x[,c(2)]) %>%
  mutate(PC3 = pcafeat$x[,c(3)]) %>%
  mutate(cluster = clusters) %>%
  mutate(major_cluster = ifelse(cluster %in% (
    as.data.frame(table(rmap_clust$cluster)) %>%
      filter(Freq >= 10) %>%
      pull(Var1) %>%
      as.numeric()
  ), cluster, 0)) %>%
  # mutate(major_cluster = factor(major_cluster,
  #                               levels = c(1, 10, 28, 30, 0),
  #                               labels = c(
  #                                 "MC-1", 
  #                                 "MC-2",
  #                                 "MC-3",
  #                                 "MC-4",
  #                                 "Other"
  #                               ))) %>%
  inner_join(rmap_samples, by = "clean_name")
rmap_clust <- inner_join(rmap_clust, rmap_clustfull[,c(which(colnames(rmap_clustfull) == "clean_name"), 
                                                       which(! colnames(rmap_clustfull) %in% colnames(rmap_clust)))])

gg <- fviz_pca_var(pcafeat, col.var="black", alpha.var = .5, 
                   select.var = list(cos2 = .15), 
                   repel = TRUE, scale. = 4.5) +
  geom_point(data = rmap_clust, 
             mapping = aes(x = PC1, y = PC2, color = Condition),
             alpha = .5, size = 3.5) +
  # scale_y_continuous(limits = c(-6, 4.5)) +
  # scale_x_continuous(limits = c(-6, 6)) +
  theme_bw(base_size = 18)
gg$layers <- rev(gg$layers)
gg$labels$x <- gsub(gg$labels$x, pattern = "Dim1", replacement = "PC1")
gg$labels$y <- gsub(gg$labels$y, pattern = "Dim2", replacement = "PC2")
biplt <- gg + labs(title = "PCA Biplot (+ Hierarchical clustering)", color = "Cluster") 
biplt + ggsave(filename = "results/biplot_clust.png", height = 7, width = 11)

gg <- fviz_pca_var(pcafeat, col.var="black", alpha.var = .5, 
                   select.var = list(cos2 = .15), 
                   repel = TRUE, scale. = 4.5) +
  geom_point(data = rmap_clust, 
             mapping = aes(x = PC1, y = PC2, color = as.factor(major_cluster)),
             alpha = .5, size = 3.5) +
  # scale_y_continuous(limits = c(-6, 4.5)) +
  # scale_x_continuous(limits = c(-6, 6)) +
  theme_bw(base_size = 18)
gg$layers <- rev(gg$layers)
gg$labels$x <- gsub(gg$labels$x, pattern = "Dim1", replacement = "PC1")
gg$labels$y <- gsub(gg$labels$y, pattern = "Dim2", replacement = "PC2")
biplt <- gg + labs(title = "PCA Biplot (+ Hierarchical clustering)", color = "Major Cluster") 
biplt + ggsave(filename = "results/biplot_majorclust.png", height = 7, width = 11)
library(viridis)

# Exon
gg <- fviz_pca_var(pcafeat, col.var="black", alpha.var = .5, 
                   select.var = list(cos2 = .15), 
                   repel = TRUE, scale. = 4.5) +
  geom_point(data = rmap_clust, 
             mapping = aes(x = PC1, y = PC2,
                           color = as.numeric(as.character(`Exon__LogP enrichment (+values depleted)`))),
             alpha = .5, size = 3.5) +
  scale_y_continuous(limits = c(-6, 4.5)) +
  scale_x_continuous(limits = c(-6, 6)) +
  theme_bw(base_size = 18)
gg$layers <- rev(gg$layers)
gg$labels$x <- gsub(gg$labels$x, pattern = "Dim1", replacement = "PC1")
gg$labels$y <- gsub(gg$labels$y, pattern = "Dim2", replacement = "PC2")
biplt <- gg + labs(title = "PCA Biplot (+ Hierarchical clustering)",
                   color = "Exon \nLog2(Obs/Exp)") +
  scale_color_viridis(begin = 0, end = .95, option = "C")
biplt
biplt + ggsave(filename = "results/biplot_exon.png", height = 7, width = 11)




library(ggsci)

condition_cmap <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "RdBu"))(12)[c(1:4, 9:12)]
names(condition_cmap) <- c("RNH", "IgG", "ACTD", "WKKD", "delta-HC", "FLAG", 'D210N', "S9.6")
condpca <- rmap_clust %>%
  arrange(desc(match(Condition, names(condition_cmap)))) %>%
  ggplot(mapping = aes(x = PC1, y = PC2, 
                       color = factor(Condition,
                                      levels = names(condition_cmap)), 
                       text = sample_name)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_point(size = 3, alpha = 1) +
  scale_y_continuous(limits = c(-6, 4.5)) +
  scale_x_continuous(limits = c(-6, 6)) +
  theme_bw(base_size = 18) +
  theme(title = element_text(size = 20)) +
  xlab(biplt$labels$x) +
  ylab(biplt$labels$y) +
  scale_color_manual(values = condition_cmap) +
  labs(title = "Condition", color = "Condition") 
condpca + ggsave(filename = "results/condition_full_pca.png", 
                 height = 6.5, width = 11)


condpca <- rmap_clust %>%
  arrange(desc(match(Condition, names(condition_cmap)))) %>%
  ggplot(mapping = aes(x = PC1, y = PC2, 
                       color = factor(Condition,
                                      levels = names(condition_cmap)), 
                       text = sample_name)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_point(size = 3, alpha = 1) +
  scale_y_continuous(limits = c(-6, 6)) +
  scale_x_continuous(limits = c(-8, 6)) +
  theme_bw(base_size = 18) +
  theme(title = element_text(size = 20)) +
  xlab(biplt$labels$x) +
  ylab(biplt$labels$y) +
  scale_color_manual(values = condition_cmap) +
  labs(title = "Condition", color = "Condition") +
  geom_text_repel(mapping = aes(x = PC1, y = PC2, 
                                label = ifelse(PC1 < 0 & Condition %in% c("RNH", "IgG", 
                                                                          "ACTD", "WKKD"),
                                               clean_name, "")), 
                  show.legend = FALSE, color = "black", force = 13)
condpca + ggsave(filename = "results/condition_full_labeled_pca.png", 
                 height = 6.5, width = 11)


modepca <- rmap_clust %>%
  ggplot(mapping = aes(x = PC1, y = PC2, 
                       color = as.factor(mode), 
                       text = sample_name)) +
  geom_point(size = 3.5) +
  theme_bw(base_size = 15) +
  xlab(paste0("PC1 (", pct_pc1, "% of variance)")) +
  ylab(paste0("PC2 (", pct_pc2, "% of variance)")) +
  scale_color_discrete(name = "Mode") +
  labs(title = "DRIP Type") 









# library(plotly)
# 
# fig <- plot_ly(rmap_clust, x = ~PC1, y = ~PC2, z = ~PC3, text = ~sample_name, color = ~as.factor(cluster))
# fig <- fig %>% add_markers()
# fig

uncertain <- c("SRX6427715" = 0, 
               "SRX6427722" = 0,
               "SRX5849996" = 1,
               "SRX7583979" = 0,
               "SRX6686232" = 0,
               "SRX6686234" = 0,
               "SRX6686235" = 0,
               "SRX6686233" = 0,
               "SRX1025907" = 0,
               "SRX5696400" = 0,
               "SRX5696401" = 0,
               "SRX5696402" = 0,
               "SRX5696403" = 0,
               "SRX5696404" = 0,
               "SRX5696405" = 0)
new_pass <- names(uncertain)[which(uncertain == 1)]


rmap_clust2 <- rmap_clust %>%
  mutate(pass = ifelse(cluster == 1 | SRX %in% new_pass, "pass", "fail"))
table(rmap_clust2$pass, rmap_clust2$Condition)
plt <- rmap_clust2 %>%
  ggplot(mapping = aes(x = PC1, y = PC2, 
                       color = pass, 
                       text = sample_name)) +
  geom_point() 
plt
plt + ggsave(filename = "results/label_pca.png", height = 5, width = 8)


clust_lab <- ifelse(rmap_clust2$pass == "pass", 1, 0)
smp_size <- floor(.75 * length(clust_lab))
set.seed(42)
train_ind <- sample(seq_len(length(clust_lab)), size = smp_size)

train_x <- rmap_features[train_ind,]
train_y <- clust_lab[train_ind]
test_x <- rmap_features[-train_ind,]
test_y <- clust_lab[-train_ind]

## XGBoost ##
library(xgboost)
# Check with validation set
dtrain <- xgb.DMatrix(data = train_x, label = train_y)
cv <- xgb.cv(data = dtrain, nrounds = 6, nthread = 2, nfold = 5,
             metrics = list('error'),
             max_depth = 5, eta = 1, objective = "binary:logistic")
cv
# Good -- train model on train data
bst <- xgboost(data = train_x, label = train_y, max.depth = 5,
               eta = 1, nthread = 2, nrounds = 6, objective = "binary:logistic")
bst
# Test on test data
pred <- predict(bst, test_x)
prediction <- as.numeric(pred > .80)
which(test_y != prediction)
names(pred) <- names(test_y)
pred <- pred[order(pred)]
pred


# 100% correct predictions. Good model. Show tree + importance
importance_matrix <- xgb.importance(model = bst)
print(importance_matrix)
png("results/xgb_importance_barchart.png", height = 5, width = 7, units = "in", res = 400)
par(mar = c(3, 12, 3, 3))
xgb.plot.importance(importance_matrix = importance_matrix, left_margin = NULL)
dev.off()

# Rebuild based on these features
xgbfeatures <- importance_matrix$Feature
# Check with validation set
dtrain <- xgb.DMatrix(data = train_x[, xgbfeatures], label = train_y)
cv <- xgb.cv(data = dtrain, nrounds = 6, nthread = 2, nfold = 5,
             metrics = list('error'), 
             max_depth = 5, eta = 1, objective = "binary:logistic")
cv
# Good -- train model on train data
bst <- xgboost(data = train_x[, xgbfeatures], label = train_y, max.depth = 5,
               eta = 1, nthread = 2, nrounds = 6, objective = "binary:logistic")
bst
# Test on test data
pred <- predict(bst, test_x[, xgbfeatures])
prediction <- as.numeric(pred > .80)
which(test_y != prediction)
names(pred) <- names(test_y)
pred <- pred[order(pred)]
pred

# Save the model
xgb.save(bst, fname = "analysis/xgb_DRIP_group_HS_binary_11_17_2020.model")
xgb.save(bst, fname = "../RSeq/helpers/data/xgb_DRIP_group_HS_binary_11_17_2020.model")
# Save features
save(xgbfeatures, file = "analysis/xgb_DRIP_group_HS_binary_11_17_2020.features.rda")
save(xgbfeatures, file = "../RSeq/helpers/data/xgb_DRIP_group_HS_binary_11_17_2020.features.rda")

## Linear Regression for scoring ##
train_x2 <- as.data.frame(train_x)
colnames(train_x2) <- gsub(colnames(train_x2), pattern = " |\\(|\\)|\\/", replacement = "_")
colnames(train_x2) <- gsub(colnames(train_x2), pattern = "^([0-9]+)", replacement = "d\\1")
test_x2 <- as.data.frame(test_x)
colnames(test_x2) <- gsub(colnames(test_x2), pattern = " |\\(|\\)|\\/", replacement = "_")
colnames(test_x2) <- gsub(colnames(test_x2), pattern = "^([0-9]+)", replacement = "d\\1")
library(caret)
train.control <- trainControl(method = "cv", number = 5)
# Train the model -- good R2
train_x2$train_y <- train_y
model <- train(form = train_y ~ 
                 corr_median + 
                 Intron__Log2_Ratio__obs_exp_ +
                 SINE__Log2_Ratio__obs_exp_ + 
                 LINE__Log2_Ratio__obs_exp_ + 
                 pct_aligned, 
               data = train_x2, method = "lm",
               trControl = train.control)
model
summary(model)
model$resample
# Test model -- good predictions 
lmpred <- predict(model, newdata = test_x2)
round(lmpred, 2)
lmpred <- round(lmpred[order(lmpred)] * 100, 2)
lmpred

# Get formula
cc <- model$finalModel$coefficients
(eqn <- paste("Y =", paste(round(cc[1],3), paste(round(cc[-1],3), names(cc[-1]), sep=" * ", collapse=" + "), sep=" + "), "+ e"))
plot(model$finalModel)

# Save model 
lmfit <- model$finalModel
save(lmfit, file = "analysis/lm_DRIP_group_HS_nonbinary_11_17_2020.rda")
save(lmfit, file = "../RSeq/helpers/data/lm_DRIP_group_HS_nonbinary_11_17_2020.rda")

# Get the full dataset predictions
rmap_features2 <- rmap_features
colnames(rmap_features2) <- gsub(colnames(rmap_features2), pattern = " |\\(|\\)|\\/", replacement = "_")
colnames(rmap_features2) <- gsub(colnames(rmap_features2), pattern = "^([0-9]+)", replacement = "d\\1")
lmfull <- predict(model, newdata = rmap_features2)

pred <- predict(bst, rmap_features[, xgbfeatures])
names(pred) <- rownames(rmap_features)
pred2 <- pred[order(pred)]
pred2

clust_lab2 <- clust_lab[order(pred)]
prediction <- as.numeric(pred2 > 0.80)
which(clust_lab2 != prediction) # Only one wrong!
names(pred) <- names(clust_lab)

rmap_quality <- data.frame(
  clean_name = names(clust_lab),
  cluster = clust_lab,
  xgbprob = pred,
  xgbverdict = as.numeric(pred > 0.80),
  lmscore = lmfull * 100
)

write_csv(rmap_quality, path = "../RSeq/helpers/data/rmap_quality_DRIP_11_17_2020.csv")


# Save all features & quality
rmap_scored_full <- rmap_quality %>%
  right_join(y = rmap_full, by = "clean_name")
write_csv(rmap_scored_full, path = "results/scored_rmap_full_11_17_2020.csv")


# Plot relationship between XGB and LM
library(ggrepel)
ggplot(rmap_quality, aes(x = lmscore, y = xgbprob, 
                         color = factor(xgbverdict,
                                        labels = c("FAIL", "PASS")))) +
  geom_point(size = 3.5) +
  ggrepel::geom_label_repel(mapping = aes(x = lmscore, y = xgbprob, 
                                          label = ifelse(xgbprob > .2 & xgbprob < .8, 
                                                         as.character(clean_name), "")),
                            color = "dimgrey", size = 4, nudge_y = -.05,
                            segment.alpha = .5, box.padding = .5, force = 3) +
  scale_x_continuous(limits = c(-50, 150)) +
  labs(title = "Linear regression vs XGB classifier") +
  ylab("XGB probability of PASS") +
  xlab("Linear model sample score") +
  geom_hline(aes(yintercept = .8), color='black', linetype='dashed', alpha = .5) +
  scale_color_manual(values = c("firebrick", "forestgreen"), name = "Verdict") +
  annotate("text", x = -30, y = .8, label = "PASS/FAIL boundary", color = "black", vjust = -.5) +
  theme_bw(base_size = 15) +
  ggsave(filename = "results/lm_vs_xgb.png", height = 6, width = 9.6)













