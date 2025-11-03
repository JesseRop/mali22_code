#!/usr/bin/env nextflow

// Wrapper script to run get_estimates_from_umi_counts.py with YASCP default parameters
// Based on the yascp cellbender module defaults

params.input_tenx = '/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/{cellranger_runs,cellranger_runs_rmHsapiens}/Pf_all_genes/5736STDY11771536/outs/raw_feature_bc_matrix' // change as needed - path to 10x directories
params.script = '/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code_lnk/cell_calling/get_estimates_from_umi_counts.py'
params.id_decode = '/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/pf_solo_mixed_mdata_decode.csv'


// Default parameters from YASCP cellbender module
params.estimate_params_umis = [
    expected_nemptydroplets_umi_cutoff: 0,
    method_estimate_ncells: "dropletutils::barcoderanks::inflection",  
    lower_bound_umis_estimate_ncells: 1000,
    method_estimate_nemptydroplets: "dropletutils::barcoderanks::knee",
    lower_bound_umis_estimate_nemptydroplets: 10,
    upper_bound_umis_estimate_nemptydroplets: 250,
    estimate_nemptydroplets_umi_add_factor: 0,
    estimate_nemptydroplets_subtract_cell_factor: 0,
    min_droplets: 0
]


// Create decode channel from CSV
id_decode_ch = Channel
    .fromPath(params.id_decode)
    .splitCsv(header: true)
    .map { row -> tuple(row.irods_id, row.sample_nm) }

// Channel: find 10x matrix directories
mtx_dirs = Channel
    .fromPath(params.input_tenx, type: 'dir', checkIfExists: true)
    .map { dir -> 
        def parts = dir.toString().split('/')
        def run_type = parts.find { it == 'cellranger_runs' || it == 'cellranger_runs_rmHsapiens' }
        def study_id = parts[parts.findIndexOf { it.contains('5736STDY') }]
        tuple(study_id, dir, run_type)
    } // Extract STDY ID and run type
    .combine(id_decode_ch, by: 0) // Join with decode channel on STDY ID
    .map { study_id, dir, run_type, msc_id -> tuple(msc_id, dir, run_type) } // Remap to use MSC ID

mtx_dirs.view()

// Process: run UMI estimates for each 10x directory
process GET_UMI_ESTIMATES {
    tag { "${sample_id} - ${run_type}" }
    cpus 1
    memory '20 GB'
    time '1.h'

    publishDir { 
        def base = '/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/'
        run_type == 'cellranger_runs' ? "${base}/${sample_id}/cbender_custom_wHsPf/umi_estimates" : "${base}/${sample_id}/cbender_custom_rmHs/umi_estimates"
    },
    mode: 'copy',
    pattern: "umi_count_estimates*",
    overwrite: true

    input:
        tuple val(sample_id), path(mtx_dir), val(run_type)

    output:
        tuple val(sample_id), 
            val(run_type),
            path("*.txt"),
            path("*.tsv.gz"),
            path("*.png")

    script:
    """
    # Activate the virtualenv for this job
    source /software/team222/jr35/yascp_utils/yascp_ncells_estimate/bin/activate

    python ${params.script} \\
        --tenxdata_path ${mtx_dir} \\
        --output_file umi_count_estimates \\
        --expected_nemptydroplets_umi_cutoff ${params.estimate_params_umis.expected_nemptydroplets_umi_cutoff} \\
        --method_estimate_ncells ${params.estimate_params_umis.method_estimate_ncells} \\
        --lower_bound_umis_estimate_ncells ${params.estimate_params_umis.lower_bound_umis_estimate_ncells} \\
        --method_estimate_nemptydroplets ${params.estimate_params_umis.method_estimate_nemptydroplets} \\
        --lower_bound_umis_estimate_nemptydroplets ${params.estimate_params_umis.lower_bound_umis_estimate_nemptydroplets} \\
        --upper_bound_umis_estimate_nemptydroplets ${params.estimate_params_umis.upper_bound_umis_estimate_nemptydroplets} \\
        --estimate_nemptydroplets_add_umifactor ${params.estimate_params_umis.estimate_nemptydroplets_umi_add_factor} \\
        --estimate_nemptydroplets_subtract_dropletfactor ${params.estimate_params_umis.estimate_nemptydroplets_subtract_cell_factor} \\
        --estimate_nemptydroplets_min_nemptydroplets ${params.estimate_params_umis.min_droplets}
    """
}

workflow {
    // Run UMI estimation on all input directories
    umi_estimates = GET_UMI_ESTIMATES(mtx_dirs)
    
    // Log completion with run type information
    umi_estimates.subscribe { 
        result ->
        def (sample_id, run_type, txt_files, gz_files, png_files) = result
        log.info """
        [Summary] Completed sample: ${sample_id}
                 Run type: ${run_type}
                 Output location: ${run_type == 'cellranger_runs' ? 'cbender_custom_wHsPf' : 'cbender_custom_rmHs'}
                 Generated files: ${txt_files.size()} txt, ${gz_files.size()} gz, ${png_files.size()} png
        """
    }
}