#!/usr/bin/env nextflow

// module load nextflow-23.10.0 

// params.id_decode = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/irods_id_sk21_mdata_cln.csv"
params.id_decode = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/pf_solo_mixed_mdata_decode.csv"  

// - Standard scripts for running soupx
params.scrpt = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/seurat_sce/soupx.R"

// - non corrected objects
params.input_seu = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/seu_obj/bpc_seu.RDS"

// - raw objects
// params.input_raw_mtx = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf/5736STDY*/outs/raw_feature_bc_matrix/"
params.input_raw_mtx = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf/5736STDY*/outs/raw_feature_bc_matrix.h5"

// - souporcell output folder
params.hvg = 500

// - souporcell output folder
// params.o_file_name = "soupx_filt_seu_norm.RDS"
params.o_sufx = "soupx_sct"

// - compute resources for first process
ncores="4"
mem="20 GB"


// - Create channels 
input_ch = Channel
				.fromPath(params.input_seu)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file, file.getParent(), file.baseName.replaceAll(/.RDS/, "") + "_" + "${params.o_sufx}" + ".RDS") }

// input_ch.view()

input_raw_ch = Channel
				.fromPath(params.input_raw_mtx)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file) }
                // .groupTuple( by: 0 ) necessary when using Read10X

// input_raw_ch.view()

// - Create channels 
id_ch = Channel
				.fromPath(params.id_decode)
				.splitCsv( header: true )
				.map { file -> tuple( file.irods_id, file.sample_nm ) }

// id_ch.view()

input_all_ch = id_ch
                    .combine(input_raw_ch, by:0)
                    .map { file -> tuple(file[1], file[0], file[2]) }
                    .combine(input_ch, by:0)
                    .map { file -> tuple(file[0], file[3], file[2], file[4], file[5]) }

input_all_ch.view()

// - Run soupx
process SOUPX {
    memory "${mem}"
    cpus "${ncores}"

    errorStrategy 'ignore' 

    tag "Soupx on ${sample_nm} reads"

    // publishDir params.outdir_index_temp, mode: "copy"
    publishDir "${out_p}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), path(seu), path(raw_mtx), val(out_p), val(o_file_name)
	
    output:
    tuple val(sample_nm), path("*")

	script:
	"""
    module load HGI/softpack/groups/team222/Pf_scRNAseq/16

    Rscript ${params.scrpt} ${seu} ${raw_mtx} ${o_file_name} ${params.hvg} 

    
	"""			
    
}


workflow {
    
    soupx_out_ch = SOUPX(input_all_ch)
    soupx_out_ch.view()


}