# Setup
library(methylKit)
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)
library(SummarizedExperiment)
library(GenomicRanges)

# Get the coverage directory path from the command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Error: Must provide paths for coverage directory and covariates file.", call. = FALSE)
}
coverage_dir <- args[1]
covars_file <- args[2] # New argument for covars file

# Read covariates using the argument. Nextflow stages this file to the workdir.
covars <- read.csv(covars_file, stringsAsFactors = FALSE)

covars$sample_id <- as.character(covars$sample_id)

head(covars)

# Make sure the directory exists
if (!dir.exists(coverage_dir)) {
  stop(paste("Error: Coverage directory not found:", coverage_dir), call. = FALSE)
}

# Define path to coverage files
#file.list <- list.files(coverage_dir, pattern = ".txt$", full.names = TRUE)

# Read methylation data using methRead
#myobjDB <- methRead(
#  location = as.list(file.list),
#  sample.id = as.list(as.character(covars$sample_id)),
#  assembly = "hg38",
#  treatment = covars$treatment,
#  context = "CpG",
#  dbtype = "tabix",
#  dbdir = "methylDB",
#  pipeline = "bismarkCytosineReport",
#  header = FALSE
#)


file.list <- file.path(coverage_dir, paste0(covars$sample_id, "_pe.CpG_report.txt"))
# Constructs file paths explicitly in the order of covars$sample_id

# Check if all files exist (critical step!)
if (!all(file.exists(file.list))) {
  missing <- file.list[!file.exists(file.list)]
  stop("These files are missing:\n", paste(missing, collapse = "\n"))
}

file.list.txt<- unlist(file.list)  # flatten list to character vector

#  Make sure your file list is a list of strings (1 string per element)
print(file.list.txt)


# If file.list.txt is a nested list, flatten it first
file.list.fixed <- as.list(unlist(file.list.txt))  # ensures list of 120 single string elements

length(file.list.fixed)  # should be 120


# Make sure sample.id is a list of strings of the same length
sample.id.fixed <- as.list(as.character(covars$sample_id))  # convert to character then list

length(sample.id.fixed)  # should be 120

#  Make sure treatment is a simple integer or numeric vector matching length 120
treatment.fixed <- covars$treatment  # should already be int vector of length 120

length(treatment.fixed)  # confirm 120

# 5. Read methylation data into methylKit
myobjDB <- methRead(
  location = file.list.fixed,
  sample.id = sample.id.fixed,
  assembly = "hg38",
  treatment = treatment.fixed,
  context = "CpG",
  dbtype = "tabix",
  dbdir = "methylDB",
  pipeline = "bismarkCytosineReport",
  header = FALSE
)



stopifnot(all.equal(getSampleID(myobjDB), covars$sample_id))

# Optional checks
print(getSampleID(myobjDB))         # Check methylKit object sample names
print(covars$sample_id)           # Should match exactly


# Sample IDs in methylKit object
getSampleID(myobjDB)

# Should match your covariates
all.equal(getSampleID(myobjDB), covars$sample_id)  # should return TRUE


# Filter samples by coverage
filtered.myobj <- filterByCoverage(myobjDB, lo.count = 10, lo.perc = NULL,
                                 hi.count = NULL, hi.perc = 99.9)

# Unite filtered samples
meth.min <- methylKit::unite(filtered.myobj, min.per.group = 10L)


# plots
pdf("clusterSamples_meth.min.pdf", h=10, w=16) 
clusterSamples(meth.min, dist="correlation", method="ward", plot=TRUE) 
dev.off()

pdf("PCAplot_meth.min.pdf", h=10, w=10) 
PCASamples(meth.min, screeplot=TRUE) 
dev.off()

pdf("PCA_meth.min.pdf", h=10, w=10)
PCASamples(meth.min)
dev.off()

# Save the unified object and covariates for the next step
save(myobjDB, meth.min, covars, filtered.myobj, file = "data_prepped.RData")
