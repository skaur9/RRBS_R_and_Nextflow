# This script performs DML analysis with adjustment for sex and saves output to a specific directory.

# Setup
library(methylKit)

# Set the output directory
output_dir <- "results/dml_sex_adjust"

# Load the prepped data
load("data_prepped.RData")

# Ensure covariates match sample order in meth.min
design.df <- covars[match(getSampleID(meth.min), covars$sample_id), ]

# Model: treatment + sex
myDiff_sex <- calculateDiffMeth(
  meth.min,
  covariates = data.frame(sex = design.df$sex),
  overdispersion = "MN",
  test = "Chisq",
  mc.cores = 8
)

# Get significant DMLs
sigDMLs_sex <- getMethylDiff(myDiff_sex, difference = 10, qvalue = 0.01)

# Write results to file in the output directory
write.table(getData(sigDMLs_sex), file = file.path(output_dir, "differential_methylation_results_sex_adjust.tsv"), sep = "\t", quote = FALSE)

# Save the object
save(myDiff_sex, file = file.path(output_dir, "dml_sex_adjust_results.RData"))
