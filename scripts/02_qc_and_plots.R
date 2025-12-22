# Setup
library(methylKit)
library(ggplot2)
library(corrplot)
library(pheatmap)
library(GenomicRanges)

# Load prepped data
load("data_prepped.RData")

# Ensure covariates match sample order in meth.min
design.df <- covars[match(getSampleID(meth.min), covars$sample_id), ]

## Plots
# Cluster samples
pdf("clusterSamples_meth.min.pdf", h = 10, w = 16)
clusterSamples(meth.min, dist = "correlation", method = "ward", plot = TRUE)
dev.off()


# PCA analysis (using methylKit's function)
pdf("PCAplot_meth.min.pdf", h = 5, w = 5)
PCASamples(meth.min, screeplot = TRUE)
dev.off()

pdf("PCA_meth.min.pdf", h = 10, w = 10)
PCASamples(meth.min)
dev.off()

## Covariates vs PercentMethylation:
meth.mat <- percMethylation(meth.min)
non_zero_var <- apply(meth.mat, 1, function(x) var(x, na.rm = TRUE) > 0)
meth.mat.filtered.clean <- meth.mat[non_zero_var, ]
meth.mat.filtered.clean <- meth.mat.filtered.clean[complete.cases(meth.mat.filtered.clean), ]
pca <- prcomp(t(meth.mat.filtered.clean), scale. = TRUE)

pca.df <- as.data.frame(pca$x)
pca.df$sex <- as.factor(design.df$sex)
pca.df$BMI_mother <- design.df$BMI_mother
pca.df$maternal_age <- design.df$mage
pca.df$BMI_child <- design.df$BMI_child
pca.df$gest_age <- design.df$ga
pca.df$treatment <- as.factor(design.df$treatment)

pdf("PCA_covariates.pdf", h = 10, w = 16)
ggplot(pca.df, aes(PC1, PC2, color = treatment)) + geom_point() + theme_minimal() + ggtitle("PCA of Methylation ~ treatment")
ggplot(pca.df, aes(PC1, PC2, color = sex)) + geom_point() + theme_minimal() + ggtitle("PCA of Methylation ~ sex")
ggplot(pca.df, aes(PC1, PC2, color = gest_age)) + geom_point() + theme_minimal() + ggtitle("PCA of Methylation ~ gest_age")
ggplot(pca.df, aes(PC1, PC2, color = maternal_age)) + geom_point() + theme_minimal() + ggtitle("PCA of Methylation ~ maternal_age")
ggplot(pca.df, aes(PC1, PC2, color = BMI_mother)) + geom_point() + theme_minimal() + ggtitle("PCA of Methylation ~ BMI_mother")
ggplot(pca.df, aes(PC1, PC2, color = BMI_child)) + geom_point() + theme_minimal() + ggtitle("PCA of Methylation ~ BMI_child")
dev.off()

# Correlation plot of covariates

# Create a correlation matrix to identify highly correlated variables and exclude them from downstream analysis

cor.mat <- cor(design.df[, c("sex", "mage", "BMI_mother", "BMI_child", "ga", "alcohol", "diab_mother", "diab_child", "rygning1", "rygning2")], use = "complete.obs")
corrplot::corrplot(cor.mat, method = "circle")


# Plot All numeric variables to idendtify corrlated variables 

# 1. Select numeric covariates starting from column 9
numeric.covs <- design.df[, 9:ncol(design.df)]

# Drop variables with more than 50% missing
numeric.covs <- numeric.covs[, colMeans(is.na(numeric.covs)) < 0.5]

# Drop rows with too many missing values
numeric.covs <- numeric.covs[rowMeans(is.na(numeric.covs)) < 0.5, ]

# 4. Compute correlation matrix
cor.mat <- cor(numeric.covs, use = "pairwise.complete.obs")

# 5. Function to compute p-values
cor.p <- function(mat) {
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  for (i in 1:n) {
    for (j in 1:n) {
      test <- suppressWarnings(cor.test(mat[, i], mat[, j]))
      p.mat[i, j] <- test$p.value
    }
  }
  p.mat
}
p.mat <- cor.p(numeric.covs)

pdf("corrplot.pdf", h = 10, w = 16)
# Plot
corrplot::corrplot(cor.mat, method = "circle")


# 6. draw the upper triangle of the correlation matrix (to avoid duplicating mirrored values).
# Color correlations by strength (blue = negative, red = positive).
# hide insignificant correlations

corrplot::corrplot(
  cor.mat,
  method = "color",
  type = "upper",
  order = "hclust",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 0.8,
  number.cex = 0.6,
  p.mat = p.mat,        # matrix of p-values
  sig.level = 0.05,     # significance cutoff
  insig = "blank",        # mark insignificant as blank
  pch.cex = 2,          # size of the significance symbol
  pch.col = "black",    # color of the significance symbol
  col = colorRampPalette(c("blue", "white", "red"))(200)
)
dev.off()

# Save the current state 
save.image(file="qc_and_plots.RData")
