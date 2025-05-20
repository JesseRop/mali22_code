#!/usr/bin/env nextflow

// module load nextflow-23.10.0 

// NOTE - running scdblfinder for a variety of settings

// - non corrected objects
params.input_seu = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/seu_obj/bpc_seu.RDS"

// - output directory
// params.o_dir= "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/seu_obj/"

// - souporcell output folder
params.mthd = "HarmonyIntegration"
params.mthd_nm = "harmony"
// params.umap_name = "UMAP_"
params.cluster_resoln = 1
// params.cl_nm = "harmony_clusters"
// params.umap_redctn_nm = "umap.harmony"
// params.new_redctn = "integrated.harmony"
params.intgrtn_pc = 30
params.integrtd_pc_man = "yes"
// params.elbo_nm = "Elbow_pc_plot"
params.findpc_mthd = "perpendicular line"
params.min_dist = 0.3
params.nneighbs = 30


// - compute resources for first process
params.ncores="3"
params.mem="40 GB"


// - Standard scripts for running souporcell
params.scrpt = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/seurat_sce/seu_integrtn_pc_auto.R"

// - Create channels 

input_ch = Channel
				.fromPath(params.input_seu)
				// .map { file -> tuple(file.getParent().toString().split('\\/')[13], file, file.getParent()) }
                .map { file -> tuple(file.getParent().toString().split('\\/')[13], file, file.getParent(), file.baseName+ "_" + "${params.mthd_nm}" + ".RDS", file.baseName+ "_" + "${params.mthd_nm}" + "_elbo.RDS") }


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
    // errorStrategy 'ignore'

    // errorStrategy 'retry' // MSC66 giving error

    tag "seurat normalization and umap on ${sample_nm} reads"

    // publishDir params.outdir_index_temp, mode: "copy"
    publishDir "${out_p}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), path(seu), val(out_p), val(nm), val(nm_elbo)
	
    output:
    tuple val(sample_nm), path(nm), path(nm_elbo)

	script:
	"""
    module load HGI/softpack/groups/team222/Pf_scRNAseq/19

    Rscript ${params.scrpt} ${seu} ${params.mthd} ${nm} ${params.cluster_resoln} "${params.mthd_nm}_clusters" "umap.${params.mthd_nm}" "integrated.${params.mthd_nm}" ${params.intgrtn_pc} ${params.integrtd_pc_man} ${nm_elbo} "${params.findpc_mthd}" ${params.min_dist} ${params.nneighbs}
    
	"""			
    
}


workflow {
    
    seu_out_ch = SEU_NORMALIZATION(input_ch)
    seu_out_ch.view()



}