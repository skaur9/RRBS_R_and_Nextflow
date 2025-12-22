# This script performs DML analysis with adjustment for sex and saves output to a specific directory.

# Setup
library(methylKit)

# 1. Read the new, explicit environment variable
num.cores <- as.numeric(Sys.getenv("R_NUM_CORES")) 

# 2. Add a sanity check to prevent the NA error
if (is.na(num.cores) || num.cores < 1) {
    num.cores <- 1 # Default to 1 if not set
}

# 3. Register the cores
library(BiocParallel)
register(MulticoreParam(workers = num.cores))

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
  test = "F",
  suffix = "sex_adjust"
  #mc.cores = 12
)

# Get significant DMLs
sigDMLs_sex <- getMethylDiff(myDiff_sex, difference = 10, qvalue = 0.01)

myDiff25p_sigDMLs_sex <- getMethylDiff(sigDMLs_sex, difference = 25, qvalue = 0.01)

# Write results to file in the output directory
write.table(getData(sigDMLs_sex), file = file.path(output_dir, "differential_methylation_results_sex_adjust.tsv"), sep = "\t", quote = FALSE)

write.table(getData(myDiff25p_sigDMLs_sex), file = file.path(output_dir, "myDiff25p_sigDMLs_sex.txt"), sep = "\t", quote = FALSE)

# Visualize the distribution of hypo/hyper-methylated bases/regions per chromosome.
pdf(file.path(output_dir, "diffMethPerChr_myDiff_sex.pdf"), h = 10, w = 16)
diffMethPerChr(myDiff_sex, plot = TRUE, qvalue.cutoff = 0.01, meth.cutoff = 25, exclude = c("GL000194.1", "GL000195.1", "GL000220.1", "GL000225.1", "MT", "Y"), keep.empty.chrom = FALSE)
dev.off()

# Save the object
save(myDiff_sex, file = file.path(output_dir, "dml_sex_adjust_results.RData"))
