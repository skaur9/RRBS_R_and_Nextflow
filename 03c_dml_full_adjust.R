# This script performs DML analysis with a full set of covariates and saves output to a specific directory.

# Setup
library(methylKit)

# Set the output directory
output_dir <- "results/dml_full_adjust"

# Load the prepped data
load("data_prepped.RData")

# Ensure covariates match sample order in meth.min
design.df <- covars[match(getSampleID(meth.min), covars$sample_id), ]

# Model: treatment + sex + mage + BMI_mother + ga + BMI_child
myDiff <- calculateDiffMeth(
  meth.min,
  covariates = data.frame(
    sex = design.df$sex,
    maternal_age = design.df$mage,
    BMI_mother = design.df$BMI_mother,
    gest_age = design.df$ga,
    BMI_child = design.df$BMI_child
  ),
  overdispersion = "MN",
  test = "Chisq",
  mc.cores = 8
)

# Get significant DMLs
sigDMLs <- getMethylDiff(myDiff, difference = 10, qvalue = 0.01)
myDiff25p <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01)

# Write results to files in the output directory
write.table(getData(sigDMLs), file = file.path(output_dir, "differential_methylation_results_full_adjust.tsv"), sep = "\t", quote = FALSE)
write.table(getData(myDiff25p), file = file.path(output_dir, "myDiff25p.txt"), sep = "\t", quote = FALSE)

# Visualize the distribution of hypo/hyper-methylated bases/regions per chromosome.
pdf(file.path(output_dir, "diffMethPerChr.pdf"), h = 10, w = 16)
diffMethPerChr(myDiff, plot = TRUE, qvalue.cutoff = 0.01, meth.cutoff = 25, exclude = c("GL000194.1", "GL000195.1", "GL000220.1", "GL000225.1", "MT", "Y"), keep.empty.chrom = FALSE)
dev.off()

# Save the object
save(myDiff, sigDMLs, file = file.path(output_dir, "dml_full_adjust_results.RData"))
