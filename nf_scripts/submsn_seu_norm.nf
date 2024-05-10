#!/usr/bin/env nextflow

// module load nextflow-23.10.0 

// NOTE - running scdblfinder for a variety of settings

// - non corrected objects
params.input_seu = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/seu_obj/bpc_seu.RDS"

// - output directory
// params.o_dir= "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/seu_obj/"

// - souporcell output folder
params.o_file_name = "seu_norm.RDS"
params.umap_name = "UMAP_"
params.hvg = 500
params.cluster_resoln = 0.4

// - compute resources for first process
ncores="4"
mem="100 GB"


// - Standard scripts for running souporcell
params.scrpt = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/seurat_sce/seu_stdnorm.R"

// - Create channels 

input_ch = Channel
				.fromPath(params.input_seu)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file, file.getParent()) }

// o_dir_ch = Channel
// 				.fromPath(params.o_dir)
// 				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file) }

// input_op_ch = input_ch
//                     .combine(o_dir_ch, by:0)

input_ch.view()

// - Run souporcell assuming K=1 - remapping, variant calling and vartrix done only once and reused for other Ks
process SEU_NORMALIZATION {
    memory "${mem}"
    cpus "${ncores}"

    // errorStrategy 'ignore' // MSC66 giving error

    tag "doublet finder on ${sample_nm} reads"

    // publishDir params.outdir_index_temp, mode: "copy"
    publishDir "${out_p}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), path(seu), val(out_p)
	
    output:
    tuple val(sample_nm), path("${params.o_file_name}")

	script:
	"""
    module load HGI/softpack/groups/team222/Pf_scRNAseq/11

    Rscript ${params.scrpt} ${seu} ${params.hvg} ${params.o_file_name} ${params.umap_name} ${params.cluster_resoln}

    
	"""			
    
}


workflow {
    
    seu_out_ch = SEU_NORMALIZATION(input_ch)
    seu_out_ch.view()

    // bcodes_bam_soupc_ch = bcodes_bam_al_ch
    //     .map { file -> tuple(file[0], file[1], file[2], file[3], file[5], file[6]) }
    //     .combine(soupc1_ch, by:[0,1])
    // bcodes_bam_soupc_ch.view()
    
    // soupc2plus_ch = SOUPC_2PLUS(soupc1_ch, k_ch)
    // // soupc2plus_ch = SOUPC_2PLUS(bcodes_bam_soupc_ch, k_ch)
    // // soupc2plus_ch.view()
    // soupc2plus_ll_ch = soupc2plus_ch
    //     .map { file -> tuple(file[0], file[1], file[5]) }
    //     .groupTuple(by: [0,1])
    
    // soupc2plus_ll_ch.view()

    // LL_KNEE_PLOT(soupc2plus_ll_ch)

}