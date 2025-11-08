#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { UMI_ESTIMATES } from './cbender_/umi_estimates'
include { CELLBENDER } from './cbender_/cellbender'
include { METADATA_CSV } from './cbender_/metadata_csv'

// Default parameters
params.outdir = '/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data'
params.rawdir = '/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data'
params.input_tenx = '/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/{cellranger_runs,cellranger_runs_rmHsapiens}/Pf_all_genes/5736STDY*/outs/raw_feature_bc_matrix'
params.decode_csv = '/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/pf_solo_mixed_mdata_decode.csv'
params.test_sample = null  // Add test sample parameter
    
// Scripts and environments
params.script_umi = '/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code_lnk/cell_calling/get_estimates_from_umi_counts.py'
params.cellbender_script = '/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts_lnk/cbender_cellNo_edropsNo.sh'
params.metadata_script = '/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code_lnk/python_nbooks/cbender_mdata_csv.py'
params.yascp_utils = '/software/team222/jr35/yascp_utils'
params.cellbender_env = '/software/team222/jr35/cellbender_gpu/cellbender_gpu_env'
    
// UMI estimation parameters
params.estimate_params_umis = [
    expected_nemptydroplets_umi_cutoff: 0,
    method_estimate_ncells: "dropletutils::barcoderanks::inflection",  
    lower_bound_umis_estimate_ncells: 251,
    method_estimate_nemptydroplets: "dropletutils::barcoderanks::inflection,dropletutils::barcoderanks::knee,0.33",
    lower_bound_umis_estimate_nemptydroplets: 10,
    upper_bound_umis_estimate_nemptydroplets: 250,
    estimate_nemptydroplets_umi_add_factor: 0,
    estimate_nemptydroplets_subtract_cell_factor: 0,
    min_droplets: 0
]

// CellBender parameters
// params.l_rates = ['0.0001', '0.000005']
// params.epochs = ['150', '250']

params.l_rates = "0.0001,0.000005"
params.epochs = "150,250"

l_rates_list = (params.l_rates instanceof List) ? params.l_rates : params.l_rates.toString().split(',')
epochs_list = (params.epochs instanceof List) ? params.epochs : params.epochs.toString().split(',')

// l_rates_ch = Channel.from( l_rates_list )
// epochs_ch = Channel.from( epochs_list )


params.ncores = 4
params.mem = 60000
// params.species_cranger = ['cellranger_runs_rmHsapiens', 'cellranger_runs']
params.species_cranger = ['cellranger_runs_rmHsapiens']
// params.species_cbender = ['_rmHs', '_wHsPf']
params.species_cbender = ['_rmHs']


// Metadata CSV parameters
params.obs_fname = 'cb_mdata.csv'
params.var_fname = 'cb_gene_mdata.csv'


// Create channels
id_decode_ch = Channel
    .fromPath(params.decode_csv)
    .splitCsv(header: true)
    .filter { row -> 
        // If test_sample is specified, only process that sample
        params.test_sample ? row.irods_id == params.test_sample : true 
    }
    .map { row -> tuple(row.irods_id, row.sample_nm) }

// id_decode_ch.view

mtx_dirs = Channel
    .fromPath(params.input_tenx, type: 'dir', checkIfExists: true)
    .map { dir -> 
        def parts = dir.toString().split('/')
        def run_type = parts.find { it == 'cellranger_runs' || it == 'cellranger_runs_rmHsapiens' }
        def irods_id = parts[parts.findIndexOf { it.contains('5736STDY') }]
        tuple(irods_id, dir, run_type)
    }
    .combine(id_decode_ch, by: 0)
    .map { irods_id, dir, run_type, msc_id -> tuple(irods_id, dir, run_type, msc_id) }

// mtx_dirs.view()

workflow {
    println "l_rates: ${l_rates_list}"
    println "epochs: ${epochs_list}"
    
    // Run UMI estimates
    UMI_ESTIMATES(mtx_dirs)
    // UMI_ESTIMATES.out.umi_results.view()
    
    // Prepare input for CellBender
    cellbender_input = UMI_ESTIMATES.out.umi_results
        .flatMap { irods_id, run_type, exp_cells, exp_drops, gz_files, png_files, msc_id ->
            def cb_suffix = run_type == 'cellranger_runs' ? '_wHsPf' : '_rmHs'
            def raw_mtx = "${params.rawdir}/${run_type}/Pf_all_genes/${irods_id}/outs/raw_feature_bc_matrix.h5"
            def combinations = []
            
            l_rates_list.each { lr ->
                epochs_list.each { ep ->
                    // Normalize lr to a plain decimal string without trailing zeros (e.g. '0.0001')
                    // def lr_raw = lr.toString()
                    // def lr_str = lr_raw.contains('.') ? new BigDecimal(lr_raw).stripTrailingZeros().toPlainString() : lr_raw
                    // Epochs are integers (e.g. '150', '250') — keep as-is
                    // def ep_str = ep.toString()

                    def o_dir = "${params.outdir}/processed/Pf/${msc_id}/cbender_custom2${cb_suffix}/cb_LR${lr}_E${ep}/"
                    def o_file = "cb.h5"
                    combinations << [msc_id, irods_id, cb_suffix, raw_mtx, o_dir, o_file, lr, ep, exp_cells, exp_drops]
                }
            }
            return combinations
        }
            
    
    cellbender_input.view()

    // Run CellBender
    CELLBENDER(cellbender_input)
    
    // Prepare input for metadata CSV generation
    metadata_input = CELLBENDER.out.cb_outputs
        .map { file -> 
            tuple(file[0], file[1], file[2], file[3])
        }
    
    metadata_input.view()

    // Generate metadata CSVs
    mdata_ch = METADATA_CSV(metadata_input)

    mdata_ch.view()
}
