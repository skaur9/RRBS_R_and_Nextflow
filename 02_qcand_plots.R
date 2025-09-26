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
numeric.covs <- design.df[, c("mage", "BMI_mother", "BMI_child", "ga")]
cor.mat <- cor(numeric.covs, use = "pairwise.complete.obs")
p.mat <- cor.p(numeric.covs) # You'll need to define the cor.p function
corrplot::corrplot(cor.mat, method = "circle")
corrplot::corrplot(cor.mat, method = "color", type = "upper", order = "hclust", addCoef.col = "black", tl.col = "black", p.mat = p.mat, sig.level = 0.05, insig = "blank")

# Save the current state if needed
save.image(file="qc_and_plots.RData")
