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
params.input_raw_mtx = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf_all_genes/5736STDY*/outs/raw_feature_bc_matrix.h5"

// - output directory
params.o_path = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/"

// - souporcell output folder
params.o_dir = "edrops"
params.o_prfx = "wo_retain"

params.limt = ["50", "100", "150", "200"]
limt_v = (params.limt instanceof List) ? params.limt : params.limt.toString().split(',')
limt_v_bcranks = ['1'] + limt_v.flatten()

limt_ch = Channel.from(limt_v)
limt_bcranks_ch = Channel.from(limt_v_bcranks)

// limt_ch.view()
// limt_bcranks_ch.view()

params.retain_nm = ["wo_retain", "w_retain"]
params.retain_val = ["Inf", "NULL"]

nms = (params.retain_nm instanceof List) ? params.retain_nm : params.retain_nm.toString().split(',')
vals = (params.retain_val instanceof List) ? params.retain_val : params.retain_val.toString().split(',')


retain_ch = Channel.from( [nms, vals].transpose() )

// retain_ch.view()

// - fdr threshold to consider cell as significantly a cell
params.fdr_l= "0.01"

// - compute resources for first process
ncores="2"
mem="10 GB"


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


input_all_ch = id_ch
                    .combine(input_raw_ch, by:0)
                    .combine(limt_ch) 
                    .combine(retain_ch)
                    .map { file -> tuple(file[1], file[2], "${params.o_path}"+file[1]+"/edrops", file[3], file[5], "edrops_"+file[3]+"_"+file[4]+".RDS", "edrops_"+file[3]+"_"+file[4]+"_slim.RDS") }


// input_all_ch.view()

input_bcrank_ch = id_ch
                    .combine(input_raw_ch, by:0)
                    .combine(limt_bcranks_ch) 
                    .map { file -> tuple(file[1], file[2], "${params.o_path}"+file[1]+"/edrops", file[3], "bcrank_"+file[3]+".RDS", "bcrank_"+file[3]+"_knee.RDS", "bcrank_"+file[3]+"_inflctn.RDS", "bcrank_"+file[3]+"_uniq.RDS") }

// input_bcrank_ch.view()

qc_stats_ch = id_ch
                    .combine(input_raw_ch, by:0)
                    .map { file -> tuple(file[1], file[2], "${params.o_path}"+file[1]+"/edrops", "qc_stats.RDS") }

qc_stats_ch.view()


// - Run edrops
process EDROPS {
    memory "${mem}"
    cpus "${ncores}"

    errorStrategy 'ignore' 

    tag "edrops on ${sample_nm} reads"

    // publishDir params.outdir_index_temp, mode: "copy"
    publishDir "${out_p}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), path(raw_p), val(out_p), val(limt), val(retain), val(out_nm), val(out_nm_slim)
	
    output:
    tuple val(sample_nm), path("*")

	script:
	"""
    module load HGI/softpack/groups/team222/Pf_scRNAseq/32

    Rscript ${params.scrpt} ${raw_p} ${limt} ${retain} ${params.fdr_l} ${out_nm}  ${out_nm_slim}
    
	"""			
    
}

ncores_bc="1"
mem_bc="10 GB"

// - Run edrops
process BCRANK {
    memory "${mem_bc}"
    cpus "${ncores_bc}"

    // errorStrategy 'ignore' 

    tag "bcranks on ${sample_nm} reads"

    // publishDir params.outdir_index_temp, mode: "copy"
    publishDir "${out_p}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), path(raw_p), val(out_p), val(limt), val(out_bcrank), val(out_bcrank_knee), val(out_bcrank_inflctn), val(out_bcrank_uniq)
	
    output:
    tuple val(sample_nm), path("*")

	script:
	"""
    module load HGI/softpack/groups/team222/Pf_scRNAseq/32

    Rscript ${params.bcrank_scrpt} ${raw_p} ${limt} ${out_bcrank} ${out_bcrank_knee} ${out_bcrank_inflctn} ${out_bcrank_uniq}
    
	"""			
    
}

// - Run edrops
process QC_STATS {
    memory "${mem_bc}"
    cpus "${ncores_bc}"

    // errorStrategy 'ignore' 

    tag "qc stats on ${sample_nm} reads"

    // publishDir params.outdir_index_temp, mode: "copy"
    publishDir "${out_p}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), path(raw_p), val(out_p), val(out_qc_stats)
	
    output:
    tuple val(sample_nm), path("*")

	script:
	"""
    module load HGI/softpack/groups/team222/Pf_scRNAseq/32

    Rscript ${params.qc_stats_scrpt} ${raw_p} ${out_qc_stats}
    
	"""			
    
}

workflow {
    
    edrops_out_ch = EDROPS(input_all_ch)
    // edrops_out_ch.view()

    bcrank_out_ch = BCRANK(input_bcrank_ch)
    // bcrank_out_ch.view()

    stats_out_ch = QC_STATS(qc_stats_ch)
    // stats_out_ch.view()

}


