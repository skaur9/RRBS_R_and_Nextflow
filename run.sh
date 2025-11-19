#!/bin/bash

module load nextflow/24.04.5
module load singularity

nextflow run main.nf \
    --outdir results \
    -c config/pbspro.config,config/custom.config \
    -resume
