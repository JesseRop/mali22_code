#!/usr/bin/env nextflow

// module load nextflow-23.10.0 

// NOTE - running scdblfinder for a variety of settings

// - non corrected objects
params.filt_obj = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC33/seu_obj/filt_seu_norm.RDS"

// - soupx corrected objects
params.spx_obj = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC33/seu_obj/filt_seu_norm_supx.RDS" 

// - output name
// params.o_file = "scdblf_outs.RDS"
params.o_sufx = "scdblf"


// - compute resources for first process
ncores="4"
mem="40 GB"

// - Standard scripts for running souporcell
params.scrpt = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/seurat_sce/scdblfindr.R"


// - Create channels 

filt_ch = Channel
				.fromPath(params.filt_obj)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file.baseName, file) }

spx_ch = Channel
				.fromPath(params.spx_obj)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file.baseName.replaceAll(/_supx/, ""), file, file.getParent(), file.baseName.replaceAll(/.RDS/, "") + "_" + "${params.o_sufx}" + ".RDS") }
				// .map { file -> tuple( "${file.irods_id}, ${file.sample_nm}" ) } NB - this does not work

// filt_ch.view()
// spx_ch.view()

all_objs_ch = filt_ch
                    .combine(spx_ch, by:[0,1])

all_objs_ch.view()


// - Run souporcell assuming K=1 - remapping, variant calling and vartrix done only once and reused for other Ks
process SCDBLFINDER_RUNS {
    memory "${mem}"
    cpus "${ncores}"

    // errorStrategy 'ignore' // MSC66 giving error

    tag "sc doublet finder on ${sample_nm} reads"

    publishDir "${out_p}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), val(obj_nm), path(filt), path(spx), val(out_p), val(o_file_name)
	
    output:
    tuple val(sample_nm), path("${o_file_name}")

	script:
	"""
    module load HGI/softpack/groups/team222/Pf_scRNAseq/11

    Rscript ${params.scrpt} ${filt} ${spx} ${o_file_name}

    
	"""			
    
}


workflow {
    
    sdblf_ch = SCDBLFINDER_RUNS(all_objs_ch)
    sdblf_ch.view()


}