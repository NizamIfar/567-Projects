#install.packages("h2o")
#install.packages("geometry")
library(dplyr)    # for data manipulation
library(ggplot2)  # for data visualization
library(tidyr)    # for data reshaping
library(geometry)
library(h2o)  # for fitting GLRMs
library(cluster)
h2o.init()

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

#normalize datasets
for(i in 1:ncol(b2)) {
  b2[, i] <- normalizeb(b2[, i], b[, i])
}
b <- apply(b, 2, normalize)

#convert datasets to H2O format  
x_train <- as.h2o(b)
x_train2 <- as.h2o(b2)

#fit autoencoder model
fit <- h2o.deeplearning(
  x = names(x_train),
  training_frame = x_train, 
  activation = "Tanh",
  autoencoder = TRUE,
  hidden = c(400, 5),  
  epochs = 1000,
  export_weights_and_biases = TRUE
)

#get layer 2 deep features
z2 <- as.matrix(h2o.deepfeatures(fit, x_train, layer = 2)) 
z2b <- as.matrix(h2o.deepfeatures(fit, x_train2, layer = 2))

# STEP 1: Select similar subset using autoencoder 
# --- after we’ve computed z2 (clinical) and z2b (real‐world) ---

# 1. Build the convex hull on the clinical deep features
u <- convhulln(z2)

# 2. Test which real‐world points lie inside it
in_hull <- inhulln(u, z2b)   # logical vector, TRUE if z2b[i,] is inside hull

# 3. Subset your real‐world data accordingly
similar_points   <- which(in_hull)
dissimilar_points <- which(!in_hull)

similar_data   <- b2[similar_points, ]
dissimilar_data <- b2[dissimilar_points, ]

cat("Points inside convex hull:   ", length(similar_points), "\n")
cat("Points outside convex hull: ", length(dissimilar_points), "\n")

# STEP 2: Cluster dissimilar data
#get the real-world data not in the hull (dissimilar data)
dissimilar_realworld <- b2[!in_hull, ]
dissimilar_z2b <- z2b[!in_hull, ]
d <- dist(dissimilar_z2b)

#check for NAs in the data
any(is.na(dissimilar_z2b))
#[1] FALSE

#using k-means clustering with optimal k determined by the elbow method
set.seed(123)  # For reproducibility

#function to calculate total within-cluster sum of squares
calculate_wss <- function(data, max_k = 20) {
  wss <- numeric(max_k)
  for (i in 1:max_k) {
    kmeans_temp <- kmeans(data, centers = i, nstart = 25, iter.max = 50)
    wss[i] <- kmeans_temp$tot.withinss
  }
  return(wss)
}

#calculate WSS for different k values
wss <- calculate_wss(dissimilar_z2b, max_k = 20)

#plot elbow method
elbow_plot <- ggplot(data.frame(k = 1:20, wss = wss), aes(x = k, y = wss)) +
  geom_line() +
  geom_point(size = 3) +
  labs(title = "Elbow Method for Optimal k",
       x = "Number of clusters (k)",
       y = "Total within-cluster sum of squares") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

print(elbow_plot)

#based on the elbow plot, choose optimal k (let's say k = 3)
k <- 3  #this can be adjusted based on the elbow plot

#perform k-means clustering with the chosen k
set.seed(123)
kmeans_result <- kmeans(dissimilar_z2b, centers = k, nstart = 25, iter.max = 100)
clusters <- kmeans_result$cluster
avg_sil_km <- mean(silhouette(kmeans_result$cluster, d)[, 3])

# 1. Prepare data.frames for plotting
cluster_df <- data.frame(
  LD1     = dissimilar_z2b[, 1],
  LD2     = dissimilar_z2b[, 2],
  cluster = factor(clusters),
  stringsAsFactors = FALSE
)

centers_df <- data.frame(
  LD1     = tapply(dissimilar_z2b[, 1], clusters, mean),
  LD2     = tapply(dissimilar_z2b[, 2], clusters, mean),
  cluster = factor(1:k),
  stringsAsFactors = FALSE
)

# 2. Plot clusters + centroids
ggplot(cluster_df, aes(x = LD1, y = LD2, color = cluster)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_point(data = centers_df, aes(x = LD1, y = LD2),
             shape = 8, size = 4, stroke = 1.2) +
  labs(title = "K-means Clusters in 2D Latent Space",
       x = "Latent Dim 1", y = "Latent Dim 2",
       color = "Cluster") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))


# 3. Compute centroids in the ORIGINAL feature space
#    and find top 5 differentiating features per cluster
centroids_orig <- aggregate(dissimilar_realworld,
                            by = list(cluster = clusters),
                            FUN = mean)

overall_mean <- colMeans(dissimilar_realworld, na.rm = TRUE)

cat("Top 5 differentiating features per cluster:\n")
for(i in 1:k) {
  # pull out the means for cluster i as a numeric vector
  cl_means <- as.numeric(centroids_orig[i, -1])
  names(cl_means) <- colnames(dissimilar_realworld)
  
  # compute difference vs overall
  diffs <- cl_means - overall_mean
  
  # pick top 5 by absolute deviation
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

# STEP 3: Compare autoencoder to PCA
#required libraries
library(vegan)     # for procrustes
library(CCA)       # for canonical correlation
library(geometry)  # for convhulln, inhulln
library(ggplot2)
library(gridExtra) # for arranging multiple plots

# 3a) Compute PCA embeddings (5 dimensions to match autoencoder bottleneck)
L <- 5
pca_res      <- prcomp(b, center = TRUE, scale. = FALSE)           # clinical PCA
pca_clinical <- pca_res$x[, 1:L]                                   # clinical scores
pca_realworld <- scale(b2, center = TRUE, scale = FALSE) %*%
  pca_res$rotation[, 1:L]                          # project real-world

# 3b) Reconstruction error comparison
#   PCA reconstruction on clinical data
rec_pca <- pca_clinical %*% t(pca_res$rotation[, 1:L])
mse_pca <- mean((b - rec_pca)^2)
#autoencoder reconstruction error (train set)
recon_ae <- as.matrix(h2o.anomaly(fit, x_train, per_feature = FALSE))
mse_ae   <- mean(recon_ae)

cat("Reconstruction MSE (PCA, 5 PCs):", round(mse_pca, 5), "\n")
cat("Reconstruction MSE (Autoencoder):", round(mse_ae, 5), "\n\n")

# 3c) Alignment of embeddings
#procrustes
proc <- procrustes(X = pca_clinical, Y = z2[, 1:L])
cat("Procrustes RMSE:", round(proc$ss, 5), "\n")

#  Canonical Correlation Analysis
cca_res <- cc(z2[, 1:L], pca_clinical)
cat("CCA correlations:", round(cca_res$cor, 4), "\n\n")

# 3d) Convex-hull agreement in 5D
hull_ae <- convhulln(z2[, 1:L])
in_ae   <- inhulln(hull_ae, z2b[, 1:L])

hull_pc <- convhulln(pca_clinical)
in_pc   <- inhulln(hull_pc, pca_realworld)

agreement_table <- table(Autoencoder = in_ae, PCA = in_pc)
print("Convex-hull classification agreement:")
print(agreement_table)
cat("\nPercent agreement:",
    round(100 * sum(diag(agreement_table)) / sum(agreement_table), 2), "%\n\n")

# 3e) Variance explained by first 2 PCs
var_exp <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
cat("Var explained PC1:", round(var_exp[1],1), "%; ",
    "PC2:", round(var_exp[2],1), "%; ",
    "Total (PC1+PC2):", round(sum(var_exp[1:2]),1), "%\n\n")

# 3f) 2×2 Visualization: AE vs PCA, clinical vs real-world
#    Prepare data.frames
ae_clin_df <- data.frame(D1 = z2[,2], D2 = z2[,1], Source = "AE-Clinical")
pca_clin_df <- data.frame(D1 = pca_clinical[,2], D2 = pca_clinical[,1], Source = "PCA-Clinical")
ae_rw_df   <- data.frame(D1 = z2b[,2], D2 = z2b[,1], Source = "AE-Real-world")
pca_rw_df  <- data.frame(D1 = pca_realworld[,2], D2 = pca_realworld[,1], Source = "PCA-Real-world")

plot_df <- rbind(ae_clin_df, pca_clin_df, ae_rw_df, pca_rw_df)

p <- ggplot(plot_df, aes(x = D1, y = D2)) +
  geom_point(alpha = 0.4, size = 1) +
  facet_wrap(~ Source, ncol = 2) +
  labs(x = "Dimension 1", y = "Dimension 2") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"))

print(p)

# STEP 4: Evaluate model performance on real‐world data

# 0) load libraries
library(caret)    # for confusionMatrix()
library(pROC)     # for ROC/AUC

# 1) read in your real‐world labels  
real_world_labels <- readRDS("~/Downloads/Real_world_ABPL.RDS")
# ensure it's a factor with levels "0","1"
test_y <- factor(real_world_labels, levels = c(0,1))

# 2) prepare your training DF
train_df <- data.frame(b, Outcome = factor(y, levels = c(0,1)))

# 3) fit a logistic model on clinical data
glm_mod <- glm(Outcome ~ ., data = train_df, family = binomial)

# 4) prepare test DF
test_df <- data.frame(b2)  # features only

# 5) get predicted probabilities on real-world
pred_prob <- predict(glm_mod, newdata = test_df, type = "response")

# 6) ROC curve & AUC
roc_res <- roc(test_y, pred_prob)
cat("AUC:", round(auc(roc_res),3), "\n")

# 7) pick optimal threshold
opt <- coords(roc_res, "best",
              ret = c("threshold","sensitivity","specificity"),
              best.method = "youden")

thresh <- as.numeric(opt["threshold"])
sens   <- as.numeric(opt["sensitivity"])
spec   <- as.numeric(opt["specificity"])
cat("Youden threshold:", round(thresh,3),
    " sens:", round(sens,3),
    " spec:", round(spec,3), "\n")

# 8) binarize at that cutpoint
pred_class <- factor(ifelse(pred_prob >= thresh, 1, 0), levels = c(0,1))

# 9) final confusion matrix
cm <- confusionMatrix(pred_class, test_y, positive = "1")
print(cm)

# Optionally, 10) view ROC
plot(roc_res, main="ROC on Real-World Data")
abline(v=thresh, col="red", lty=2)

# STEP 5: Interpretation and visualization 

# 0) Load required libraries
library(caret)       # train(), confusionMatrix()
library(ranger)      # fast RF backend for caret
library(pROC)        # ROC, AUC, coords()
library(xgboost)     # gradient boosting
library(randomForest) # for direct RF (if desired)
set.seed(123)

# 1) Read in your real‐world outcome labels
real_raw <- readRDS("~/Downloads/Real_world_ABPL.RDS")
test_y   <- factor(real_raw,
                   levels = c(0,1),
                   labels = c("neg","pos"))

# 2) Build your clinical‐study training set
#    (assumes 'b' is your clinical features matrix, 'y' your clinical binary outcome)
train_df <- data.frame(
  b,
  Outcome = factor(y,
                   levels = c(0,1),
                   labels = c("neg","pos"))
)

# 3) Prepare a test features data.frame
test_df <- data.frame(b2)   # your real-world features

# ----------------------------
# 4A) Logistic Regression (GLM)
# ----------------------------
glm_mod <- glm(Outcome ~ ., data = train_df, family = binomial)

# Predict probabilities on real-world
pred_prob_glm <- predict(glm_mod, newdata = test_df, type = "response")

# ROC / AUC
roc_glm <- roc(test_y, pred_prob_glm)
cat("GLM AUC:", round(auc(roc_glm), 3), "\n")

# Youden’s J cutpoint
opt_glm <- coords(roc_glm, "best",
                  ret         = c("threshold","sensitivity","specificity"),
                  best.method = "youden")
th_glm   <- as.numeric(opt_glm["threshold"])
sen_glm  <- as.numeric(opt_glm["sensitivity"])
spe_glm  <- as.numeric(opt_glm["specificity"])
cat(sprintf("GLM Youden threshold = %.3f (sens = %.3f, spec = %.3f)\n",
            th_glm, sen_glm, spe_glm))

# Confusion matrix at Youden
pred_class_glm <- factor(ifelse(pred_prob_glm >= th_glm, "pos", "neg"),
                         levels = c("neg","pos"))
cm_glm <- confusionMatrix(pred_class_glm, test_y, positive = "pos")
print(cm_glm)

# Make sure your labels are a factor with levels c(0,1)
test_labels <- factor(real_world_labels, levels = c(0,1))

# Define the four subsets by index
subsets <- list(
  Similar  = similar_points,
  Cluster1 = which(!in_hull & clusters == 1),
  Cluster2 = which(!in_hull & clusters == 2),
  Cluster3 = which(!in_hull & clusters == 3)
)

# Function to evaluate one subset
evaluate_subset <- function(name, idx) {
  Xsub   <- as.data.frame(b2[idx, , drop = FALSE])
  ysub   <- test_labels[idx]
  probs  <- predict(glm_mod, newdata = Xsub, type = "response")
  rocobj <- roc(ysub, probs)
  
  # Youden’s J
  opt    <- coords(rocobj, "best",
                   ret = c("threshold","sensitivity","specificity"),
                   best.method = "youden")
  thr    <- as.numeric(opt["threshold"])
  sens   <- as.numeric(opt["sensitivity"])
  spec   <- as.numeric(opt["specificity"])
  
  # Classify at that threshold
  preds  <- factor(ifelse(probs >= thr, 1, 0), levels = c(0,1))
  cm     <- confusionMatrix(preds, ysub, positive = "1")
  
  # Print summary
  cat("====", name, "====\n")
  cat("  n =", length(idx), "  Positive rate =", round(mean(as.numeric(as.character(ysub))),3), "\n")
  cat("  AUC:", round(auc(rocobj),3), 
      "  Youden thr:", round(thr,3), 
      "  (sens", round(sens,3), "spec", round(spec,3), ")\n")
  print(cm$table)
  cat("\n")
  
  invisible(list(roc = rocobj, cm = cm))
}

# Run it for each subset
results <- lapply(names(subsets), function(nm){
  evaluate_subset(nm, subsets[[nm]])
})

# Name the results list so names(results) matches the 4 elements
names(results) <- names(subsets)

# Now build the AUC data frame
auc_df <- data.frame(
  subset = names(results),
  AUC    = sapply(results, function(res) as.numeric(auc(res$roc)))
)

# Quick check
print(auc_df)

# --- Bar‐plot of AUC by subset ---
library(ggplot2)

ggplot(auc_df, aes(x = subset, y = AUC, fill = subset)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", AUC)), 
            vjust = -0.5, size = 3.5) +
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 11, face = "bold")
  ) +
  labs(
    title = "AUROC by Sub-population",
    x     = "",
    y     = "Area Under ROC Curve (AUC)"
  )


# --- Overlay all 4 ROC curves on one plot ---
library(pROC)

cols <- c("firebrick","steelblue","darkgreen","purple")
plot(0,0,type="n", xlim=c(1,0), ylim=c(0,1),
     xlab="1 - Specificity", ylab="Sensitivity")

i <- 1
for(nm in names(results)) {
  plot.roc(results[[nm]]$roc, add = TRUE, col = cols[i], lwd = 2)
  i <- i + 1
}
legend("bottomright", legend = names(results),
       col = cols, lwd = 2, cex = 0.8)

