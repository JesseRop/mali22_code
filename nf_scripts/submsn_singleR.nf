#!/usr/bin/env nextflow

// module load nextflow-23.10.0 

// - Standard scripts for running souporcell
params.scrpt = "/nfs/users/nfs_j/jr35/multipurpose_scripts/seurat_sce/singleR_mult.R"

// - non corrected objects
params.input_seu = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC14/seu_obj/bpseu_raw_{edrops50,cbendr}_pref_sct.RDS"

// - reference objects
params.ref_obj_p= "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/m2_all/annot_refs/msc_qcd_anotd_sce_13_14.RDS"
params.marker_obj_p= "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/m2_all/annot_refs/com_markers_localMac.RDS"


// - souporcell output folder
params.ref_stage = "stageHL"

// - suffix for output files
params.o_sufx = "pred_sR"

// - compute resources for first process
params.ncores="3"
params.mem="25 GB"


// - Create channels 

input_ch = Channel
				.fromPath(params.input_seu)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], 
                                        file, 
                                        file.getParent(), 
                                        file.baseName + "_" + "sce" + ".RDS", 
                                        // file.baseName.replace("bpseu_raw", "pred").split('_')[0..-2].join('_') + ".RDS", 
                                        file.baseName + "_" + "${params.o_sufx}" + ".RDS", 
                                        params.ref_obj_p, 
                                        params.marker_obj_p) }

input_ch.view()


// - Run souporcell assuming K=1 - remapping, variant calling and vartrix done only once and reused for other Ks
process SINGLE_R {
    memory "${params.mem}"
    cpus "${params.ncores}"

    // errorStrategy 'ignore' // MSC66 giving error

    tag "singleR stage prediction on ${sample_nm}"

    // publishDir params.outdir_index_temp, mode: "copy"
    publishDir "${out_p}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), path(query_obj_p), val(out_p), val(o_sce), val(o_pred), path(ref_obj_p), path(marker_obj_p)
    
    output:
    tuple val(sample_nm), path(o_sce), path(o_pred)

    script:
    """
    module load HGI/softpack/groups/team222/Pf_scRNAseq/32

    Rscript ${params.scrpt} ${ref_obj_p} ${marker_obj_p} ${query_obj_p} ${o_sce} ${params.ref_stage} ${o_pred}
    
    """			
    
}


workflow {
    
    seu_out_ch = SINGLE_R(input_ch)
    seu_out_ch.view()

}

