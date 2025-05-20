#!/usr/bin/env nextflow

// module load nextflow-23.10.0 

// - non corrected objects
params.input_seu = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/seu_obj/seu_norm.RDS"

// - Reference SCE object to get labels from
params.sce_ref= "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/v3labmscs.ref.sce.stgasx.RDS"

// - number of cells threshold needed to pass cos sim threshold below
params.w = 3

// - scmap cosine similarity threshold
params.cos_sim = 0.5

// - output name
params.o_sufx = "prelim_anotd"

// - column with the stage labels in the reference
params.stg_labl = "stageHL_ref"

// - compute resources for first process
ncores="5"
mem="30 GB"

// - Standard scripts for running souporcell
params.scrpt = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/seurat_sce/scmap.R"

// - Create channels 

input_ch = Channel
				.fromPath(params.input_seu)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file, file.getParent(), file.baseName + "_" + "${params.o_sufx}" + ".RDS") }


input_ch.view()

// - Run souporcell assuming K=1 - remapping, variant calling and vartrix done only once and reused for other Ks
process SCMAP_LABEL_TRANSFER {
    memory "${mem}"
    cpus "${ncores}"

    // errorStrategy 'ignore' // MSC66 giving error

    tag "scmap on ${sample_nm}"

    // publishDir params.outdir_index_temp, mode: "copy"
    publishDir "${out_p}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), path(seu), val(out_p), val(o_file)
	
    output:
    tuple val(sample_nm), path("${o_file}")

	script:
	"""
    module load HGI/softpack/groups/team222/Pf_scRNAseq/11

    Rscript ${params.scrpt} ${seu} ${params.sce_ref} ${params.w} ${params.cos_sim} ${o_file} ${params.stg_labl}

    
	"""			
    
}


workflow {
    
    seu_out_ch = SCMAP_LABEL_TRANSFER(input_ch)
    seu_out_ch.view()


}