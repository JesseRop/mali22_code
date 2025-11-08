#!/usr/bin/env nextflow

// Nextflow version of cbender_submsn_cellNo_edropsNo.sh
// Runs cellbender for each sample using expected_cells and total_droplets_included from umi_estimates.nf
// Produces cb.h5 files for downstream use in cbender_mdata_csv.nf

params.decode_csv = '/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/pf_solo_mixed_mdata_decode.csv'
params.base_dir = '/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data'
params.raw_dir = '/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data'
params.cellbender_script = '/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts_lnk/cbender_cellNo_edropsNo.sh'
params.l_rates = ['0.0001', '0.000005']
params.epochs = ['150', '250']
// params.l_rates = ['0.0001']
// params.epochs = ['150']
// params.ncores = 4
// params.mem = 80000
params.ncores = 4
params.mem = 60000

// Suffixes for Pf only and Pf+Hsapiens
params.species_cranger = ['cellranger_runs_rmHsapiens', 'cellranger_runs']
params.species_cbender = ['_rmHs', '_wHsPf']

// Read decode CSV and create sample mapping channel
sample_ch = Channel
    .fromPath(params.decode_csv)
    .splitCsv(header: true)
    .filter { row -> row.irods_id == '5736STDY11771536' }
    .map { row -> tuple(row.sample_nm, row.irods_id) }

// Expand all combinations (emit each element separately using fromList)
// Build combinations without using .cross() (avoids channel/list cross issues)
// Build combinations pairing the cranger and cbender species entries by index
// e.g. 'cellranger_runs_rmHsapiens' <-> '_rmHs' and 'cellranger_runs' <-> '_wHsPf'
combos_ch = sample_ch.flatMap { tuple_sample ->
    def (sample_name, source_irods) = tuple_sample
    def out = []
    def n = Math.min(params.species_cranger.size(), params.species_cbender.size())
    for (int i = 0; i < n; i++) {
        def cr = params.species_cranger[i]
        def cb = params.species_cbender[i]
        params.l_rates.each { lr ->
            params.epochs.each { ep ->
                out << tuple(sample_name, source_irods, cr, cb, lr, ep)
            }
        }
    }
    return out
}

sample_ch.view()
combos_ch.view()

// Find all required input files for each combo
// Note: combos_ch emits tuples: [sample_name, source_irods, cranger, cbender, l_rate, epoch]
inputs_ch = combos_ch.map { sample_name, source_irods, cranger, cbender, l_rate, epoch ->
    def raw_mtx = "${params.raw_dir}/${cranger}/Pf_all_genes/${source_irods}/outs/raw_feature_bc_matrix.h5"
    def o_dir = "${params.base_dir}/processed/Pf/${sample_name}/cbender_custom3${cbender}/cb_LR${l_rate}_E${epoch}/"
    // def o_file = "${o_dir}/cb.h5"
    def o_file = "cb.h5"
    // Read estimates produced by umi_estimates pipeline placed under the sample processed dir
    def exp_cells_file = "${params.base_dir}/processed/Pf/${sample_name}/cbender_custom${cbender}/umi_estimates/umi_count_estimates-expected_cells.txt"
    def exp_edrops_file = "${params.base_dir}/processed/Pf/${sample_name}/cbender_custom${cbender}/umi_estimates/umi_count_estimates-total_droplets_included.txt"
    tuple(sample_name, source_irods, cbender, raw_mtx, o_dir, o_file, l_rate, epoch, exp_cells_file, exp_edrops_file, params.ncores, params.mem)
}

inputs_ch.view()

// Process: run cellbender for each sample/config
process RUN_CELLBENDER {
    tag { "${sample_name} - LR${l_rate}_E${epoch} ${cbender}" }
    cpus { ncores }
    memory { mem + ' MB' }
    time '5.h'
    // Submit this process to the GPU queue
    queue 'gpu-normal'

    // Add retry strategy
    errorStrategy 'ignore'
    // maxRetries 2

    // Request GPU resources when using a cluster executor (LSF)
    // clusterOptions "-q gpu-normal -gpu 'num=1:j_exclusive=yes' -M ${mem} -R \"rusage[mem=${mem}]\" -cwd ${o_dir}"
    // clusterOptions "-q gpu-normal -gpu 'num=1:j_exclusive=yes' -M ${mem} -R \"rusage[mem=${mem}]\""
    // Use clusterOptions only for GPU flags; queue and memory are controlled by Nextflow directives
    clusterOptions "-gpu 'num=1:j_exclusive=yes'"

    // publish outputs to the input file's parent directory (dynamic)
    // mode can be 'copy', 'symlink', or 'move'
    publishDir { o_dir }, mode: 'copy', pattern: "cb*", overwrite: true

    input:
        tuple val(sample_name), val(source_irods), val(cbender), path(raw_mtx), val(o_dir), val(o_file), val(l_rate), val(epoch), path(exp_cells_file), path(exp_edrops_file), val(ncores), val(mem)

    output:
        path("cb*")

    script:
    """
    mkdir -p ${o_dir}
    
    exp_cells=\$(cat ${exp_cells_file})
    exp_edrops=\$(cat ${exp_edrops_file})

    # Run the wrapper script but force direct run mode so it executes cellbender here
    ${params.cellbender_script} ${raw_mtx} ${o_file} ${l_rate} ${epoch} \$exp_cells \$exp_edrops 

    """
}

workflow {
    // To actually run cellbender, uncomment the following two lines
    cb_h5_ch = RUN_CELLBENDER(inputs_ch)
    // cb_h5_ch.view()
}
