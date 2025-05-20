#!/usr/bin/env nextflow

// module load nextflow-23.10.0 

// - Standard scripts for running souporcell
params.scrpt = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/seurat_sce/seu_stdnorm.R"

// - non corrected objects
params.input_seu = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/seu_obj/bpc_seu.RDS"

// - output directory
// params.o_dir= "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/seu_obj/"

// - souporcell output folder
params.hvg = 500
params.normlzn = "norm"
params.umap_name = "UMAP_"
params.cluster_resoln = 0.4
params.cl_nm = "NULL"
params.umap_redctn_nm = "umap"
params.redctn = "pca"
params.integrtd_pc30 = "no"
params.findpc_mthd = "perpendicular line"
params.pca_redctn_nm = "pca"
params.min_dist = 0.3

// - compute resources for first process
params.ncores="3"
params.mem="10 GB"


// - Create channels 

input_ch = Channel
				.fromPath(params.input_seu)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file, file.getParent(), file.baseName+ "_" + "${params.normlzn}" + ".RDS", file.baseName+ "_" + "${params.normlzn}" + "_elbo.RDS") }

// o_dir_ch = Channel
// 				.fromPath(params.o_dir)
// 				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file) }

// input_op_ch = input_ch
//                     .combine(o_dir_ch, by:0)

input_ch.view()

// - Run souporcell assuming K=1 - remapping, variant calling and vartrix done only once and reused for other Ks
process SEU_NORMALIZATION {
    memory "${params.mem}"
    cpus "${params.ncores}"

    // errorStrategy 'retry' // MSC66 giving error

    tag "seurat normalization and umap on ${sample_nm} reads"

    // publishDir params.outdir_index_temp, mode: "copy"
    publishDir "${out_p}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), path(seu), val(out_p), val(nm), val(nm_elbo)
	
    output:
    // tuple val(sample_nm), path(nm + "_" + "${params.normlzn}" + ".RDS"), path(nm + "_" + "${params.normlzn}" + "_elbo.RDS")
    tuple val(sample_nm), path(nm), path(nm_elbo)

	script:
	"""
    module load HGI/softpack/groups/team222/Pf_scRNAseq/16

    Rscript ${params.scrpt} ${seu} ${params.hvg} ${nm} ${params.umap_name} ${params.cluster_resoln} ${params.cl_nm} ${params.umap_redctn_nm} ${params.redctn} ${params.integrtd_pc30} ${nm_elbo} "${params.findpc_mthd}" ${params.pca_redctn_nm} ${params.min_dist}
    
	"""			
    
}


workflow {
    
    seu_out_ch = SEU_NORMALIZATION(input_ch)
    seu_out_ch.view()

}

