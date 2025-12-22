//main.nf 

//Process 1: data_prep
process data_prep {
    input:
        path covars_file
        val coverage_dir
    output:
        path 'data_prepped.RData', emit: prepped_data
    script:
        """
        Rscript ${params.datadir}/scripts/01_data_prep.R "$coverage_dir" "${covars_file.name}"
        """
}

//Process 2: QC and plots
//process qc_and_plots {
//    input:
//        path data_prepped_file
//    output:
//       path 'qc_and_plots.RData', emit: qc_results
//    script:
//        """
//        Rscript ${params.datadir}/scripts/02_qc_and_plots.R
//        """
//}

//Process 3: DMLS without adjustment
process dml_no_covars {
    cpus 8
    input:
        path data_prepped_file
    output:
        path 'results/dml_no_covars/', emit: dml_no_covars_results
    script:
        """
        mkdir -p results/dml_no_covars
        Rscript ${params.datadir}/scripts/03a_dml_no_covars.R
        """
}

//Process 4: DMLS adjusted for sex
process dml_sex_adjust {
    cpus 8
    input:
        path data_prepped_file
    output:
        path 'results/dml_sex_adjust/', emit: dml_sex_adjust_results
    script:
        """
        mkdir -p results/dml_sex_adjust
        Rscript ${params.datadir}/scripts/03b_dml_sex_adjust.R
        """
}

//Process 5: DMLS with covariates
process dml_full_adjust {
    cpus 8
    input:
        path data_prepped_file
    output:
        path 'results/dml_full_adjust/', emit: dml_full_adjust_results
    script:
        """
        mkdir -p results/dml_full_adjust
        Rscript ${params.datadir}/scripts/03c_dml_full_adjust.R
        """
}

workflow {
    // 1. Define the input channel for the covariates file
    covars_in = file(params.covars)
    
    // 2. Run data_prep, which outputs 'prepped_data'
    data_prep_channel = data_prep(covars_in, params.coverage_dir)
    
    // 3. Pass 'prepped_data' to all downstream analysis processes
    //qc_and_plots(data_prep_channel.prepped_data)
    dml_no_covars(data_prep_channel.prepped_data)
    dml_sex_adjust(data_prep_channel.prepped_data)
    dml_full_adjust(data_prep_channel.prepped_data)
}
