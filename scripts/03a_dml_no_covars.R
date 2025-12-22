# This script performs DML analysis without any covariates and saves output to a specific directory.

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
output_dir <- "results/dml_no_covars"

# Load the prepped data
load("data_prepped.RData")

# Model: Treatment only
myDiff_DMR <- calculateDiffMeth(
  meth.min,
  overdispersion = "MN",
  test = "Chisq",
  suffix = "no_covars"
  #mc.cores = 12
)

# Get significant DMLs
sigDMLs_myDiff_DMR <- getMethylDiff(myDiff_DMR, difference = 10, qvalue = 0.01)
myDiff25p_myDiff_DMR <- getMethylDiff(sigDMLs_myDiff_DMR, difference = 25, qvalue = 0.01)

# Write results to file in the output directory
write.table(getData(sigDMLs_myDiff_DMR), file = file.path(output_dir, "differential_methylation_results_no_covars.tsv"), sep = "\t", quote = FALSE)

write.table(getData(myDiff25p_myDiff_DMR), file = file.path(output_dir, "myDiff25p_myDiff_DMR.txt"), sep = "\t", quote = FALSE)

# Visualize the distribution of hypo/hyper-methylated bases/regions per chromosome.
pdf(file.path(output_dir, "diffMethPerChr_myDiff_DMR.pdf"), h = 10, w = 16)
diffMethPerChr(myDiff_DMR, plot = TRUE, qvalue.cutoff = 0.01, meth.cutoff = 25, exclude = c("GL000194.1", "GL000195.1", "GL000220.1", "GL000225.1", "MT", "Y"), keep.empty.chrom = FALSE)
dev.off()

# Save the object for potential further use
save(myDiff_DMR, file = file.path(output_dir, "dml_no_covars_results.RData"))
