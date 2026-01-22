#!/usr/bin/env nextflow

// module load nextflow-23.10.0 

// params.id_decode = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/irods_id_sk21_mdata_cln.csv"
params.id_decode = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/pf_all_id_decode.csv"  

// - Standard scripts for running soupx
params.scrpt = "/nfs/users/nfs_j/jr35/multipurpose_scripts/seurat_sce/edrops.R"
params.bcrank_scrpt = "/nfs/users/nfs_j/jr35/multipurpose_scripts/seurat_sce/bcrank_calc.R"
params.qc_stats_scrpt = "/nfs/users/nfs_j/jr35/multipurpose_scripts/seurat_sce/qc_stats.R"

// - raw objects
// ## !! NOTE - USE raw_feature_bc_matrix.h5 and not the folder raw_feature_bc_matrix/ since nextflow can't expand to all the donor paths when using the folder after specifying 5736STDY*
params.input_raw_mtx = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/{cellranger_runs,cellranger_runs_rmHsapiens}/Pf_all_genes/5736STDY*/outs/raw_feature_bc_matrix.h5"

// - output directory
params.o_path = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/"

// params.limt = ["50", "100", "150", "200"]
params.limt = "50,100,150,200"
limt_v = (params.limt instanceof List) ? params.limt : params.limt.toString().split(',')
limt_v_bcranks = ['1'] + limt_v.flatten()

limt_ch = Channel.from(limt_v)
limt_bcranks_ch = Channel.from(limt_v_bcranks)

// limt_ch.view()
// limt_bcranks_ch.view()

// params.retain_nm = ["wo_retain", "w_retain"]
// params.retain_val = ["Inf", "NULL"]
params.retain_nm = "wo_retain,w_retain"
params.retain_val = "Inf,NULL"

nms = (params.retain_nm instanceof List) ? params.retain_nm : params.retain_nm.toString().split(',')
vals = (params.retain_val instanceof List) ? params.retain_val : params.retain_val.toString().split(',')


retain_ch = Channel.from( [nms, vals].transpose() )

// retain_ch.view()

// - fdr threshold to consider cell as significantly a cell
params.fdr_l= "0.01"

// - compute resources for first process
ncores="5"
mem="30 GB"


input_raw_ch = Channel
				.fromPath(params.input_raw_mtx)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file) }
                // .groupTuple( by: 0 ) // necessary when using Read10X

// input_raw_ch.view()

// - Create channels 
id_ch = Channel
				.fromPath(params.id_decode)
				.splitCsv( header: true )
				.map { file -> tuple( file.irods_id, file.sample_nm ) }

// id_ch.view()


// Combined input channel for all analyses
combined_input_ch = id_ch
                    .combine(input_raw_ch, by:0)
                    .map { file ->
                        def sample = file[1]
                        def raw_p = file[2]
                        def suffix = raw_p.toString().contains('cellranger_runs_rmHsapiens') ? 'edrops_rmHs' : 'edrops_wHsPf'
                        def out_p = "${params.o_path}${sample}/"+suffix
                        tuple(sample, raw_p, out_p)
                    }

combined_input_ch.view()


// - Combined process for all analyses
process EDROPS_BCRANK_QC {
    memory "${mem}"
    cpus "${ncores}"
    time '2.h'

    // errorStrategy 'ignore' 

    tag "analyses on ${sample_nm} reads"

    publishDir "${out_p}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), path(raw_p), val(out_p)
	
    output:
    tuple val(sample_nm), path("*")

	script:
	"""
    module load HGI/softpack/groups/team222/Pf_scRNAseq/34

    # Run EDROPS for all limit and retain combinations
    ${limt_v.collect { limt ->
        [nms, vals].transpose().collect { retain_pair ->
            def retain_nm = retain_pair[0]
            def retain_val = retain_pair[1]
            "Rscript ${params.scrpt} ${raw_p} ${limt} ${retain_val} ${params.fdr_l} edrops_${limt}_${retain_nm}.RDS edrops_${limt}_${retain_nm}_slim.RDS"
        }.join('\n    ')
    }.join('\n    ')}

    # Run BCRANK for all limits (including '1')
    ${limt_v_bcranks.collect { limt ->
        "Rscript ${params.bcrank_scrpt} ${raw_p} ${limt} bcrank_${limt}.RDS bcrank_${limt}_knee.RDS bcrank_${limt}_inflctn.RDS bcrank_${limt}_uniq.RDS"
    }.join('\n    ')}

    # Run QC_STATS
    Rscript ${params.qc_stats_scrpt} ${raw_p} qc_stats.RDS
    
	"""			
    
}

workflow {
    
    all_out_ch = EDROPS_BCRANK_QC(combined_input_ch)
    // all_out_ch.view()

}


