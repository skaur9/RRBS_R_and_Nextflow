# Setup
library(methylKit)
library(data.table)
library(tidyverse)
library(dplyr)

# Read covariates
covars <- read.csv("covars_all.csv", stringsAsFactors = FALSE)


# Get the coverage directory path from the command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop("Error: No coverage directory path provided.", call. = FALSE)
}
coverage_dir <- args[1]

# Make sure the directory exists
if (!dir.exists(coverage_dir)) {
  stop(paste("Error: Coverage directory not found:", coverage_dir), call. = FALSE)
}

# 2. Define path to coverage files
file.list <- list.files(coverage_dir, pattern = ".txt$", full.names = TRUE)

# Read methylation data using methRead
myobjDB <- methRead(
  location = as.list(file.list),
  sample.id = as.list(as.character(covars$sample_id)),
  assembly = "hg38",
  treatment = covars$treatment,
  context = "CpG",
  dbtype = "tabix",
  dbdir = "methylDB",
  pipeline = "bismarkCytosineReport",
  header = FALSE
)

# Filter samples by coverage
filtered.myobj <- filterByCoverage(myobjDB, lo.count = 10, lo.perc = NULL,
                                 hi.count = NULL, hi.perc = 99.9)

# Unite filtered samples
meth.min <- unite(filtered.myobj, min.per.group = 10L)

# Save the unified object and covariates for the next step
save(meth.min, covars, file = "data_prepped.RData")
