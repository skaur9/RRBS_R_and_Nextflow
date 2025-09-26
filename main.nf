params.datadir = '/home/projects/cu_10031/data/raw/DNBC_methylation_2023/2025_analysis/methylkit_analysis/'

params.coverage_dir = "/home/projects/cu_10031/data/raw/DNBC_methylation_2023/2025_analysis/child_dataset"

process data_prep {
    // Pass the parameter to the R script
    input:
        path covars_file
        val coverage_dir
    output:
        path 'data_prepped.RData'
    script:
        """
        Rscript ${params.datadir}/01_data_prep.R "$coverage_dir"
        """
}


process qc_and_plots {
    input:
        path data_prepped_file
    output:
        path 'qc_and_plots.RData'
    script:
        """
        Rscript ${params.datadir}/02_qc_and_plots.R
        """
}

process dml_no_covars {
    input:
        path data_prepped_file
    output:
        path 'results/dml_no_covars/'
    script:
        """
        mkdir -p results/dml_no_covars
        Rscript ${params.datadir}/03a_dml_no_covars.R
        """
}

// ... other DML processes
