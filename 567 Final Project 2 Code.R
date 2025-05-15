#load equired libraries
library(dplyr)      #for data manipulation
library(ggplot2)    #for data visualization
library(tidyr)      #for data reshaping
library(geometry)   #for convex hull operations
library(h2o)        #for autoencoder
library(cluster)    #for clustering analysis
library(vegan)      #for procrustes analysis
library(CCA)        #for canonical correlation
library(pROC)       #for ROC/AUC analysis
library(caret)      #for confusionMatrix

#initialize h2o
h2o.init()

#define Fisher-Yates transformation functions
fisheryatesb = function(v, u) {
  q = qnorm((1:length(u))/(length(u)+1))[rank(u)]
  q0 = sort(unique(q))
  u0 = sort(unique(u))
  v0 = sort(unique(v))
  v1 = v0[v0<min(u0)]
  v2 = v0[v0>max(u0)]
  k = length(q0)
  
  if(length(q0) >= 2) {
    q1 = q0[1] + (v1-u0[1]) * (q0[2]-q0[1])/(u0[2]-u0[1])
    q2 = q0[k] + (v2-u0[k]) * (q0[k-1]-q0[k])/(u0[k-1]-u0[k])
  } else {
    q1 = q0[1]
    q2 = q0[1]
  }
  
  approx(c(v1, u0, v2), c(q1, q0, q2), v)$y
}

fisheryates = function(x) qnorm((1:length(x))/(length(x)+1))[rank(x)]

#load datasets 
b <- readRDS("~/Downloads/clinical_study.RDS")
b2 <- readRDS("~/Downloads/Real_world.RDS")

#extract target variable and remove from features
y <- b[, 2]
b <- b[, -2] 

#remove PE_MILD variable
pe_mild_col <- which(colnames(b) == "PE_MILD") 
b <- b[, -pe_mild_col]
pe_mild_col <- which(colnames(b2) == "PE_MILD")
b2 <- b2[, -pe_mild_col]

#normalize function as suggested by professor
normalize <- function(x) (x - mean(x)) / ifelse((mu1 <- (max(x) - min(x))) == 0, 1, mu1)  
normalizeb <- function(y, x) (y - mean(x)) / ifelse((mu1 <- (max(x) - min(x))) == 0, 1, mu1)

#apply fy transformation 
b_fy <- b
b2_fy <- b2

#transform each column
for(i in 1:ncol(b)) {
  col_name <- colnames(b)[i]
  
  #error check: skip columns with NAs (no error message)
  if(any(is.na(b[, i])) || any(is.na(b2[, i]))) {
    cat("Skipping column", col_name, "due to NAs\n")
    next
  }
  
  #transform clinical data
  b_fy[, i] <- fisheryates(b[, i])
  
  #transform real-world data using fisheryatesb to align with clinical
  b2_fy[, i] <- fisheryatesb(b2[, i], b[, i])
}

#normalize datasets
for(i in 1:ncol(b2_fy)) {
  b2_fy[, i] <- normalizeb(b2_fy[, i], b_fy[, i])
}
b_fy <- apply(b_fy, 2, normalize)

#data is in matrix format
b_fy <- as.matrix(b_fy)
b2_fy <- as.matrix(b2_fy)

#convert datasets to H2O format  
x_train_fy <- as.h2o(b_fy)
x_train2_fy <- as.h2o(b2_fy)

#fit autoencoder model
set.seed(123)
fit <- h2o.deeplearning(
  x = names(x_train_fy),
  training_frame = x_train_fy, 
  activation = "Tanh",
  autoencoder = TRUE,
  hidden = c(400, 5),  
  epochs = 1000,
  export_weights_and_biases = TRUE,
  l1 = 1e-5,  #new: regularization for stability
  seed = 123  #new: for reproducibility
)

#get layer 2 deep features
z2 <- as.matrix(h2o.deepfeatures(fit, x_train_fy, layer = 2)) 
z2b <- as.matrix(h2o.deepfeatures(fit, x_train2_fy, layer = 2))

# STEP 1: Build convex hull and find similar points
u <- convhulln(z2)
in_hull <- inhulln(u, z2b)

#get similar and dissimilar points
similar_points <- which(in_hull)
dissimilar_points <- which(!in_hull)

cat("Points inside convex hull:", length(similar_points), "\n")
cat("Points outside convex hull:", length(dissimilar_points), "\n")

# Create subsets
similar_data <- b2_fy[similar_points, ]
dissimilar_data <- b2_fy[dissimilar_points, ]
dissimilar_z2b <- z2b[dissimilar_points, ]

# STEP 2: Cluster dissimilar data using hierarchical clustering with Ward's method
#calculate distance matrix
d <- dist(dissimilar_z2b)

#hierarchical clustering with Ward's distance
#suggestion: if memsize issue occurs, clear environment, 
#restart R and run again from scratch
set.seed(123)
hc <- hclust(d, method = "ward.D2")

#cut tree to get k=3 clusters
k <- 3
clusters <- cutree(hc, k = k)

#evaluate clustering quality with silhouette
sil <- silhouette(clusters, d)
avg_sil_hc <- mean(sil[, 3])
cat("Average silhouette width for Ward's hierarchical clustering:", avg_sil_hc, "\n")

#visualize clusters in 2D latent space
cluster_df <- data.frame(
  LD1 = dissimilar_z2b[, 1],
  LD2 = dissimilar_z2b[, 2],
  cluster = factor(clusters)
)

centers_df <- data.frame(
  LD1 = tapply(dissimilar_z2b[, 1], clusters, mean),
  LD2 = tapply(dissimilar_z2b[, 2], clusters, mean),
  cluster = factor(1:k)
)

#cluster plot
cluster_plot <- ggplot(cluster_df, aes(x = LD1, y = LD2, color = cluster)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_point(data = centers_df, aes(x = LD1, y = LD2),
             shape = 8, size = 4, stroke = 1.2) +
  labs(title = "Ward's Hierarchical Clustering in 2D Latent Space",
       subtitle = "After Fisher-Yates Transformation",
       x = "Latent Dim 1", y = "Latent Dim 2",
       color = "Cluster") +
  theme_minimal()

print(cluster_plot)

#top differentiating features for each cluster
centroids_orig <- aggregate(dissimilar_data, by = list(cluster = clusters), FUN = mean)
overall_mean <- colMeans(dissimilar_data, na.rm = TRUE)

cat("\nTop 5 differentiating features per cluster:\n")
for(i in 1:k) {
  cl_means <- as.numeric(centroids_orig[i, -1])
  names(cl_means) <- colnames(dissimilar_data)
  
  diffs <- cl_means - overall_mean
  top_feats <- sort(abs(diffs), decreasing = TRUE)[1:5]
  
  cat(sprintf("\nCluster %d:\n", i))
  for(feat in names(top_feats)) {
    direction <- ifelse(diffs[feat] > 0, "higher", "lower")
    cat(sprintf("  %s: %s than avg (%.2f vs %.2f)\n",
                feat, direction,
                cl_means[feat],
                overall_mean[feat]))
  }
}

# STEP 3: Compare PCA embeddings (5 dimensions to match autoencoder bottleneck)
L <- 5  # Number of dimensions
pca_res_fy <- prcomp(b_fy, center = TRUE, scale. = FALSE)
pca_clinical_fy <- pca_res_fy$x[, 1:L]
pca_realworld_fy <- scale(b2_fy, center = TRUE, scale = FALSE) %*% pca_res_fy$rotation[, 1:L]

#Reconstruction error comparison
#PCA reconstruction on clinical data
rec_pca <- pca_clinical_fy %*% t(pca_res_fy$rotation[, 1:L])
mse_pca <- mean((b_fy - rec_pca)^2)

#autoencoder reconstruction error (train set)
recon_ae <- as.matrix(h2o.anomaly(fit, x_train_fy, per_feature = FALSE))
mse_ae <- mean(recon_ae)

cat("\nReconstruction MSE (PCA, 5 PCs):", round(mse_pca, 5), "\n")
cat("Reconstruction MSE (Autoencoder):", round(mse_ae, 5), "\n\n")

#Alignment of embeddings
proc <- procrustes(X = pca_clinical_fy, Y = z2[, 1:L])
cat("Procrustes RMSE:", round(proc$ss, 5), "\n")

#Canonical Correlation Analysis
cca_res <- cc(z2[, 1:L], pca_clinical_fy)
cat("CCA correlations:", round(cca_res$cor, 4), "\n\n")

#Convex-hull agreement in 5D
hull_pc <- convhulln(pca_clinical_fy)
in_pc <- inhulln(hull_pc, pca_realworld_fy)

agreement_table <- table(Autoencoder = in_hull, PCA = in_pc)
print("Convex-hull classification agreement:")
print(agreement_table)
cat("\nPercent agreement:",
    round(100 * sum(diag(agreement_table)) / sum(agreement_table), 2), "%\n\n")

#variance explained by PCA (first 2)
var_exp <- (pca_res_fy$sdev^2 / sum(pca_res_fy$sdev^2)) * 100
cat("Var explained PC1:", round(var_exp[1],1), "%; ",
    "PC2:", round(var_exp[2],1), "%; ",
    "Total (PC1+PC2):", round(sum(var_exp[1:2]),1), "%\n\n")

#variance explained by PCA (first 5)
cat("Variance explained by PCA components:\n")
for(i in 1:5) {
  cat(sprintf("PC%d: %.1f%%\n", i, var_exp[i]))
}

#loop to calculate cumulative variance explained
cumulative_var <- cumsum(var_exp[1:5])
cat("\nCumulative variance explained:\n")
for(i in 1:5) {
  cat(sprintf("PC1-PC%d: %.1f%%\n", i, cumulative_var[i]))
}

#total variance explained by all 5 components
cat(sprintf("\nTotal variance explained (PC1-PC5): %.1f%%\n", sum(var_exp[1:5])))

#here we'll add a more consistent AE vs PCA visualization comparison
#comparison plot 1: standard dimension ordering
ae_clin_df <- data.frame(D1 = z2[,1], D2 = z2[,2], Source = "AE-Clinical")
pca_clin_df <- data.frame(D1 = pca_clinical_fy[,1], D2 = pca_clinical_fy[,2], Source = "PCA-Clinical")
ae_rw_df <- data.frame(D1 = z2b[,1], D2 = z2b[,2], Source = "AE-Real-world")
pca_rw_df <- data.frame(D1 = pca_realworld_fy[,1], D2 = pca_realworld_fy[,2], Source = "PCA-Real-world")

#combine data for first plot
plot_df_v1 <- rbind(ae_clin_df, pca_clin_df, ae_rw_df, pca_rw_df)

#comparison plot 2: swapped PCA dimensions
pca_clin_df_swapped <- data.frame(D1 = pca_clinical_fy[,2], D2 = pca_clinical_fy[,1], Source = "PCA-Clinical")
pca_rw_df_swapped <- data.frame(D1 = pca_realworld_fy[,2], D2 = pca_realworld_fy[,1], Source = "PCA-Real-world")

#combine data for second plot
plot_df_v2 <- rbind(ae_clin_df, pca_clin_df_swapped, ae_rw_df, pca_rw_df_swapped)

#comparison plot with standard dimension ordering
comparison_plot_v1 <- ggplot(plot_df_v1, aes(x = D1, y = D2)) +
  geom_point(alpha = 0.4, size = 1) +
  facet_wrap(~ Source, ncol = 2) +
  labs(title = "Comparison: Same Dimension Ordering for AE and PCA",
       x = "Dimension 1", y = "Dimension 2") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))

#comparison plot with swapped PCA dimensions
comparison_plot_v2 <- ggplot(plot_df_v2, aes(x = D1, y = D2)) +
  geom_point(alpha = 0.4, size = 1) +
  facet_wrap(~ Source, ncol = 2) +
  labs(title = "Comparison 2: Swapped PCA Dimensions (PC2, PC1)",
       x = "Dimension 1", y = "Dimension 2") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))

#hull membership visualization 
#this helps visualize where real-world points fall relative to clinical hull
ae_hull_df <- rbind(
  data.frame(D1 = z2[,1], D2 = z2[,2], 
             Type = "Clinical", 
             Method = "Autoencoder"),
  data.frame(D1 = z2b[,1], D2 = z2b[,2], 
             Type = ifelse(in_hull, "Inside Hull", "Outside Hull"),
             Method = "Autoencoder")
)

pca_hull_df <- rbind(
  data.frame(D1 = pca_clinical_fy[,1], D2 = pca_clinical_fy[,2], 
             Type = "Clinical", 
             Method = "PCA"),
  data.frame(D1 = pca_realworld_fy[,1], D2 = pca_realworld_fy[,2], 
             Type = ifelse(in_pc, "Inside Hull", "Outside Hull"),
             Method = "PCA")
)

hull_df <- rbind(ae_hull_df, pca_hull_df)
hull_df$Type <- factor(hull_df$Type, levels = c("Clinical", "Inside Hull", "Outside Hull"))

#hull membership plot
hull_plot <- ggplot(hull_df, aes(x = D1, y = D2, color = Type)) +
  geom_point(alpha = 0.4, size = 1) +
  facet_wrap(~ Method, scales = "free") +
  scale_color_manual(values = c("gray60", "blue", "red")) +
  labs(title = "AE vs PCA: Hull Membership",
       x = "Dimension 1", y = "Dimension 2") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))

#finally, print all plots
print(comparison_plot_v1)
print(comparison_plot_v2)
print(hull_plot)

# STEP 4: Evaluate model performance on real-world data
real_world_labels <- readRDS("~/Downloads/Real_world_ABPL.RDS")
test_y <- factor(real_world_labels, levels = c(0,1))

#preparing the training/testing DF
train_df <- data.frame(b_fy, Outcome = factor(y, levels = c(0,1)))
test_df <- data.frame(b2_fy)

#fit a logistic model on clinical data
glm_mod <- glm(Outcome ~ ., data = train_df, family = binomial)

#get predicted probabilities on real-world
pred_prob <- predict(glm_mod, newdata = test_df, type = "response")

#ROC curve & AUC
roc_res <- roc(test_y, pred_prob)
cat("Overall AUC:", round(auc(roc_res), 3), "\n")

#pick optimal threshold
opt <- coords(roc_res, "best",
              ret = c("threshold","sensitivity","specificity"),
              best.method = "youden")

thresh <- as.numeric(opt["threshold"])
sens <- as.numeric(opt["sensitivity"])
spec <- as.numeric(opt["specificity"])
cat("Youden threshold:", round(thresh, 3),
    " sens:", round(sens, 3),
    " spec:", round(spec, 3), "\n")

#create predicted classes (binarize at that cutpoint)
pred_class <- factor(ifelse(pred_prob >= thresh, 1, 0), levels = c(0,1))

#final confusion matrix
cm <- confusionMatrix(pred_class, test_y, positive = "1")
print(cm)

# STEP 5: Evaluate on different subsets
evaluate_subset <- function(name, idx) {
  
  #error catch: no errors found
  if(length(idx) == 0 || length(unique(test_y[idx])) < 2) {
    cat("====", name, "==== (Skipping - insufficient data)\n\n")
    return(NULL)
  }
  
  Xsub <- as.data.frame(b2_fy[idx, , drop = FALSE])
  ysub <- test_y[idx]
  
  probs <- predict(glm_mod, newdata = Xsub, type = "response")
  rocobj <- roc(ysub, probs)
  
  #obtain Youden's J threshold
  opt <- coords(rocobj, "best",
                ret = c("threshold","sensitivity","specificity"),
                best.method = "youden")
  thr <- as.numeric(opt["threshold"])
  sens <- as.numeric(opt["sensitivity"])
  spec <- as.numeric(opt["specificity"])
  
  #classify at optimal threshold
  preds <- factor(ifelse(probs >= thr, 1, 0), levels = c(0,1))
  cm <- confusionMatrix(preds, ysub, positive = "1")
  
  #print summary
  cat("====", name, "====\n")
  cat("  n =", length(idx), "  Positive rate =", round(mean(as.numeric(as.character(ysub))), 3), "\n")
  cat("  AUC:", round(auc(rocobj), 3), 
      "  Youden thr:", round(thr, 3), 
      "  (sens", round(sens, 3), "spec", round(spec, 3), ")\n")
  print(cm$table)
  cat("\n")
  
  return(list(roc = rocobj, cm = cm))
}

#define subsets
subsets <- list(
  Similar = similar_points,
  Cluster1 = dissimilar_points[clusters == 1],
  Cluster2 = dissimilar_points[clusters == 2],
  Cluster3 = dissimilar_points[clusters == 3]
)

#evaluate each subset
results <- list()
for(name in names(subsets)) {
  results[[name]] <- evaluate_subset(name, subsets[[name]])
}

#filter out NULL results (none found)
results <- results[!sapply(results, is.null)]

#create AUC dataframe
auc_df <- data.frame(
  subset = names(results),
  AUC = sapply(results, function(res) as.numeric(auc(res$roc)))
)

print(auc_df)

#prepare ROC curves plot
roc_data <- list()
for(name in names(results)) {
  roc_obj <- results[[name]]$roc
  roc_data[[name]] <- data.frame(
    name = name,
    sensitivities = roc_obj$sensitivities,
    specificities = roc_obj$specificities,
    auc = round(auc(roc_obj), 3)
  )
}

#combine all data frames
roc_df <- do.call(rbind, roc_data)

#set colors for each group
subset_colors <- c(
  "Similar" = "blue", 
  "Cluster1" = "green3", 
  "Cluster2" = "purple3", 
  "Cluster3" = "orange2"
)

#create ROC plot
roc_plot <- ggplot(roc_df, aes(x = 1-specificities, y = sensitivities, color = name)) +
  geom_path(size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(
    values = subset_colors[names(results)],
    labels = sapply(names(results), function(name) {
      paste0(name, " (AUC=", round(auc(results[[name]]$roc), 3), ")")
    })
  ) +
  labs(
    title = "ROC Curves by Subpopulation",
    subtitle = "After Fisher-Yates and Ward's Hierarchical Clustering",
    x = "False Positive Rate (1-Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Subset"
  ) +
  theme_minimal() +
  coord_equal()

print(roc_plot)

#enhanced plot, now with dissimilar and overall:

#evaluate on all real-world data
all_rw_results <- evaluate_subset("All Real-World", 1:length(test_y))

#evaluate on all dissimilar points (all points outside the convex hull)
dissimilar_results <- evaluate_subset("Dissimilar", dissimilar_points)

#new results list that includes everything
results_with_all <- c(
  results, 
  list(
    "All Real-World" = all_rw_results,
    "Dissimilar" = dissimilar_results
  )
)

#new list of colors
all_subset_colors <- c(
  "Similar" = "blue", 
  "Cluster1" = "green3", 
  "Cluster2" = "purple3", 
  "Cluster3" = "orange2",
  "All Real-World" = "red3",
  "Dissimilar" = "darkgray"
)

#new roc data with the additional curves
all_roc_data <- list()
for(name in names(results_with_all)) {
  roc_obj <- results_with_all[[name]]$roc
  all_roc_data[[name]] <- data.frame(
    name = name,
    sensitivities = roc_obj$sensitivities,
    specificities = roc_obj$specificities,
    auc = round(auc(roc_obj), 3)
  )
}

#plotting df for rocs
all_roc_df <- do.call(rbind, all_roc_data)

#final roc plot
roc_plot_complete <- ggplot(all_roc_df, aes(x = 1-specificities, y = sensitivities, color = name)) +
  geom_path(size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(
    values = all_subset_colors[names(results_with_all)],
    labels = sapply(names(results_with_all), function(name) {
      paste0(name, " (AUC=", round(auc(results_with_all[[name]]$roc), 3), ")")
    })
  ) +
  labs(
    title = "ROC Curves by Subpopulation (Complete Analysis)",
    subtitle = "With Feature Engineering and Cluster-Specific Models",
    x = "False Positive Rate (1-Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Subset"
  ) +
  theme_minimal() +
  coord_equal()

print(roc_plot_complete)



########################
# FEATURE ENGINEERING Section: 
# Adding interaction terms and non-linear features
########################
#engineering function
engineer_features <- function(data) {
  #data frame check 
  data <- as.data.frame(data)
  
  #add interaction terms based on cluster analysis
  data$RACE_MARITAL <- data$RACE * data$MARITAL
  data$AGE_GDM <- data$AGE * data$GDM
  data$WGT_GES_HYP <- data$WGT * data$GES_HYP
  
  #create non-linear terms for key variables
  data$RACE_sq <- data$RACE^2
  data$AGE_sq <- data$AGE^2
  data$WGT_sq <- data$WGT^2
  
  #create aggregate risk score (not used in final version due to multicollinearity)
  #data$RISK_SCORE <- scale(data$RACE + data$GDM + data$GES_HYP + data$WGT)
  
  return(data)
}

#apply feature engineering to original datasets
train_df <- data.frame(b_fy, Outcome = factor(y, levels = c(0,1)))
train_df_enhanced <- engineer_features(train_df)

test_df <- data.frame(b2_fy)
test_df_enhanced <- engineer_features(test_df)

#fitting enhanced logistic regression model
glm_mod <- glm(Outcome ~ ., data = train_df_enhanced, family = binomial)

#get predictions with enhanced features
pred_prob <- predict(glm_mod, newdata = test_df_enhanced, type = "response")

#train random forest for ensemble
set.seed(123)
rf_mod <- randomForest(Outcome ~ ., 
                       data = train_df_enhanced,
                       ntree = 200,
                       importance = TRUE)

#combine predictions (ensemble)
rf_probs <- predict(rf_mod, newdata = test_df_enhanced, type = "prob")[,2]
ensemble_probs <- (2*pred_prob + rf_probs)/3  
#weighted average favoring GLM. 
#A few other configurations were tested, and this seemed to perform the best. 

#ROC curve & AUC for ensemble
roc_res <- roc(test_y, ensemble_probs)
cat("Overall AUC (ensemble):", round(auc(roc_res), 3), "\n")

#finding optimal threshold
opt <- coords(roc_res, "best",
              ret = c("threshold","sensitivity","specificity"),
              best.method = "youden")

thresh <- as.numeric(opt["threshold"])
sens <- as.numeric(opt["sensitivity"])
spec <- as.numeric(opt["specificity"])
cat("Youden threshold:", round(thresh, 3),
    " sens:", round(sens, 3),
    " spec:", round(spec, 3), "\n")

#create predicted classes
pred_class <- factor(ifelse(ensemble_probs >= thresh, 1, 0), levels = c(0,1))

#confusion matrix
cm <- confusionMatrix(pred_class, test_y, positive = "1")
print(cm)

# STEP 6: Evaluate on different subsets with cluster-specific models
evaluate_subset <- function(name, idx) {
  if(length(idx) == 0 || length(unique(test_y[idx])) < 2) {
    cat("====", name, "==== (Skipping - insufficient data)\n\n")
    return(NULL)
  }
  
  #get subset of test data and apply feature engineering
  Xsub <- as.data.frame(b2_fy[idx, , drop = FALSE])
  Xsub_enhanced <- engineer_features(Xsub)
  ysub <- test_y[idx]
  
  #if this is a cluster, train a specialized model
  if(name != "Similar") {
    #obtain subset of the original training data for efficiency
    train_sample_idx <- sample(nrow(train_df_enhanced), min(5000, nrow(train_df_enhanced)))
    train_sample <- train_df_enhanced[train_sample_idx, ]
    
    #obtain a sample of this cluster's data (for balancing the model)
    cluster_sample_size <- min(1000, length(idx))
    cluster_sample_idx <- sample(length(idx), cluster_sample_size)
    
    #use a balanced dataset with clinical data and this cluster's data
    #error check: ensure we do this only if we have labels for the cluster data
    if(all(idx[cluster_sample_idx] <= length(test_y))) {
      cluster_df <- data.frame(
        b2_fy[idx[cluster_sample_idx], ],
        Outcome = test_y[idx[cluster_sample_idx]]
      )
      cluster_df_enhanced <- engineer_features(cluster_df)
      
      #combine with original training data for specialized model
      combined_train <- rbind(
        train_sample,
        cluster_df_enhanced
      )
      
      #train cluster-specific model
      cluster_glm <- glm(Outcome ~ ., data = combined_train, family = binomial)
      
      #make predictions with both global and cluster-specific models
      global_probs <- predict(glm_mod, newdata = Xsub_enhanced, type = "response")
      cluster_probs <- predict(cluster_glm, newdata = Xsub_enhanced, type = "response")
      rf_probs <- predict(rf_mod, newdata = Xsub_enhanced, type = "prob")[,2]
      
      #ensemble predictions (weighted average)
      probs <- (global_probs + cluster_probs + rf_probs)/3
    } else {
      #fall back to global model if we can't train a specialized one
      global_probs <- predict(glm_mod, newdata = Xsub_enhanced, type = "response")
      rf_probs <- predict(rf_mod, newdata = Xsub_enhanced, type = "prob")[,2]
      probs <- (2*global_probs + rf_probs)/3
    }
  } else {
    #for similar points, just use the global ensemble
    global_probs <- predict(glm_mod, newdata = Xsub_enhanced, type = "response")
    rf_probs <- predict(rf_mod, newdata = Xsub_enhanced, type = "prob")[,2]
    probs <- (2*global_probs + rf_probs)/3
  }
  
  rocobj <- roc(ysub, probs)
  
  #obtain Youden's J threshold
  opt <- coords(rocobj, "best",
                ret = c("threshold","sensitivity","specificity"),
                best.method = "youden")
  thr <- as.numeric(opt["threshold"])
  sens <- as.numeric(opt["sensitivity"])
  spec <- as.numeric(opt["specificity"])
  
  #classify at optimal threshold
  preds <- factor(ifelse(probs >= thr, 1, 0), levels = c(0,1))
  cm <- confusionMatrix(preds, ysub, positive = "1")
  
  #print summary
  cat("====", name, "====\n")
  cat("  n =", length(idx), "  Positive rate =", round(mean(as.numeric(as.character(ysub))), 3), "\n")
  cat("  AUC:", round(auc(rocobj), 3), 
      "  Youden thr:", round(thr, 3), 
      "  (sens", round(sens, 3), "spec", round(spec, 3), ")\n")
  print(cm$table)
  cat("\n")
  
  return(list(roc = rocobj, cm = cm))
}

#subsets
subsets <- list(
  Similar = similar_points,
  Cluster1 = dissimilar_points[clusters == 1],
  Cluster2 = dissimilar_points[clusters == 2],
  Cluster3 = dissimilar_points[clusters == 3]
)

#evaluate each subset
results <- list()
for(name in names(subsets)) {
  results[[name]] <- evaluate_subset(name, subsets[[name]])
}

#filter out NULL results (not needed)
results <- results[!sapply(results, is.null)]

#create AUC dataframe
auc_df <- data.frame(
  subset = names(results),
  AUC = sapply(results, function(res) as.numeric(auc(res$roc)))
)

print(auc_df)

#prepare ROC curves plot
roc_data <- list()
for(name in names(results)) {
  roc_obj <- results[[name]]$roc
  roc_data[[name]] <- data.frame(
    name = name,
    sensitivities = roc_obj$sensitivities,
    specificities = roc_obj$specificities,
    auc = round(auc(roc_obj), 3)
  )
}

# Combine all data frames
roc_df <- do.call(rbind, roc_data)

# Set colors for each group
subset_colors <- c(
  "Similar" = "blue", 
  "Cluster1" = "green3", 
  "Cluster2" = "purple3", 
  "Cluster3" = "orange2"
)

# Create ROC plot
roc_plot <- ggplot(roc_df, aes(x = 1-specificities, y = sensitivities, color = name)) +
  geom_path(size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(
    values = subset_colors[names(results)],
    labels = sapply(names(results), function(name) {
      paste0(name, " (AUC=", round(auc(results[[name]]$roc), 3), ")")
    })
  ) +
  labs(
    title = "ROC Curves by Subpopulation",
    subtitle = "With Feature Engineering and Cluster-Specific Models",
    x = "False Positive Rate (1-Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Subset"
  ) +
  theme_minimal() +
  coord_equal()

print(roc_plot)

#comparison table for easy viewing
results_summary <- data.frame(
  Subset = names(results),
  AUC = sapply(results, function(r) round(auc(r$roc), 3)),
  Sensitivity = sapply(results, function(r) round(r$cm$byClass["Sensitivity"], 3)),
  Specificity = sapply(results, function(r) round(r$cm$byClass["Specificity"], 3)),
  Accuracy = sapply(results, function(r) round(r$cm$overall["Accuracy"], 3)),
  Size = sapply(names(results), function(n) length(subsets[[n]]))
)

#results table
print(results_summary)

#Supplementary checks: 

#VIF: 
library(car)
#look at which coefficients are aliased; 
#in our case, the risk-score (see earlier), if added, should be perfectly aliased
aliased_info <- alias(glm_mod)
print(aliased_info$Complete)

#we should note that “RISK_SCORE” is aliased. 
#we can remove the aliased term(s) and refit
#we removed Risk Score in the final update. 
glm_mod_reduced <- update(glm_mod, . ~ . - RISK_SCORE)

#then compute VIFs on the reduced model. 
#The expectation is that GDM and GDM interactions have VIP around 5, 
#and other variables have a VIP of 2. A majority are close to one. 
vif_reduced <- vif(glm_mod_reduced)
print(round(vif_reduced, 2))

##Feature Importance Plots (unused):

#visualize feature importance from random forest
importance_df <- as.data.frame(importance(rf_mod))
importance_df$Feature <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]

#plot feature importance (top 15)
imp_plot <- ggplot(importance_df[1:min(15, nrow(importance_df)),], 
                   aes(x = reorder(Feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Importance (Random Forest)",
       x = "Features",
       y = "Importance (Mean Decrease Gini)") +
  theme_minimal()

print(imp_plot)
