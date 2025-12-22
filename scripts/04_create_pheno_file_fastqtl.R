# This script extracts DMRs from no covariates model and creates  phenotype file for meQTL analysis 

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
load("dml_no_covars_results.Rdata")

# Extract the DMRs from the methylDiff object
DMR_list_10 <- getMethylDiff(
    myDiff,
    qvalue = 0.01,
    difference = 10,
    type = "all" )

DMR_list_25 <- getMethylDiff(
    myDiff,
    qvalue = 0.01,
    difference = 25,
    type = "all" )

# Convert the methylDiff object to a GRanges object
DMR_GRanges_10 <- as(DMR_list_10, "GRanges")

DMR_GRanges_25 <- as(DMR_list_25, "GRanges")

# Use regionCounts to get CpG counts aggregated over the DMR regions
# myobj is our unified methylBase object (meth.min)
DMR_counts_10 <- methylKit::regionCounts(
    myobj = meth.min, 
    regions = DMR_GRanges_10
)


DMR_counts_25 <- methylKit::regionCounts(
    myobj = meth.min, 
    regions = DMR_GRanges_25
)


# Get the methylation percentage (0 to 100) for each DMR in each sample
DMR_meth_perc_10 <- percMethylation(DMR_counts_10)
DMR_meth_perc_25 <- percMethylation(DMR_counts_25)

# Note: Some tools prefer 0-1 (beta value), so we  might need to divide by 100 later.

# 1. Get the sample IDs (column names)
sample_IDs <- getSampleID(DMR_counts_25)

# 2. Get the DMR coordinates (GRanges)
DMR_coords <- getData(DMR_counts_25)[, 1:3]

# 3. Create the required 6-column prefix
# ProbeID and Group_ID can both be the Chr:Start:End coordinate for simplicity
DMR_prefix <- data.frame(
    # Column 1-3: Coordinates
    "#Chr" = DMR_coords$chr,
    "Start" = DMR_coords$start,
    "End" = DMR_coords$end,
    # Column 4: ProbeID/DMRID (Using the combined coordinate)
    "ProbeID" = paste(DMR_coords$chr, DMR_coords$start, DMR_coords$end, sep = "_"),
    # Column 5: Group_ID
    "Group_ID" = paste(DMR_coords$chr, DMR_coords$start, DMR_coords$end, sep = "_"),
    # Column 6: Strand (Use '.' as strand is often ignored for DMRs)
    "Strand" = "."
)

# 4. Combine the prefix with the methylation percentages (converting to 0-1)
FASTQTL_DMR_Matrix <- cbind(DMR_prefix, DMR_meth_perc_25 / 100)
colnames(FASTQTL_DMR_Matrix)[7:ncol(FASTQTL_DMR_Matrix)] <- sample_IDs

# 5. Save the file
write.table(
    FASTQTL_DMR_Matrix,
    file = "FASTQTL_DMR_Phenotype_25.bed",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

