process UMI_ESTIMATES {
    tag { "${irods_id}(${msc_id}) - ${run_type}" }
    cpus 1
    memory '20 GB'
    time '1.h'

    publishDir { 
        def base = "${params.outdir}/processed/Pf/"
        run_type == 'cellranger_runs' ? "${base}/${msc_id}/cbender_custom_wHsPf/umi_estimates" : "${base}/${msc_id}/cbender_custom_rmHs/umi_estimates"
    },
    mode: 'copy',
    pattern: "umi_count_estimates*",
    overwrite: true

    input:
        tuple val(irods_id), path(mtx_dir), val(run_type), val(msc_id)

    output:
        tuple val(irods_id), 
            val(run_type),
            path("umi_count_estimates-expected_cells.txt"),
            path("umi_count_estimates-total_droplets_included.txt"),
            path("*.tsv.gz"),
            path("*.png"),
            val(msc_id), emit: umi_results

    script:
    """
    source ${params.yascp_utils}/yascp_ncells_estimate/bin/activate

    python ${params.script_umi} \\
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