#!/bin/bash

# --- 1. Set Environment and Paths ---

# Define the output directory 
OUTDIR="./results" 
COVARS_FILE="/home/projects/cu_10031/data/raw/DNBC_methylation_2023/2025_analysis/methylkit_analysis/covars_all.csv"
COVERAGE_DIR="/home/projects/cu_10031/data/raw/DNBC_methylation_2023/2025_analysis/child_dataset/"

# --- 2. Load Nextflow/Container Modules ---

module load jdk/23
module load nextflow/24.10.3
module load singularity/4.3.0


# --- 3. Run Nextflow Pipeline ---

nextflow run main.nf \
    -c config/custom.config \
    -profile cluster \
    --outdir "${OUTDIR}" \
    --covars "${COVARS_FILE}" \
    --coverage_dir "${COVERAGE_DIR}" \
    -with-trace "${OUTDIR}/trace.txt" \
    -with-report "${OUTDIR}/report.html" \
    -with-timeline "${OUTDIR}/timeline.html" \
    -with-dag "${OUTDIR}/pipeline_dag.png" \
    -resume


