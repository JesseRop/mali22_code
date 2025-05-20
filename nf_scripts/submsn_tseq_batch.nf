#!/usr/bin/env nextflow

// module load nextflow-23.10.0 
params.w_dir = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/pseudotime_tradeseq_ip/sympto/m2_asex_sympto_bal"

// - non corrected objects
// params.sce_obj = "${params.w_dir}_sce.RDS"

// // - pseudotime object
// params.ptime_obj= "${params.w_dir}_ptime.RDS"

// // -  cell weights object
// params.cw_obj= "${params.w_dir}_cw.RDS"

// // -  cell weights object
// params.u_don= "${params.w_dir}_u_don.RDS"


// - number of cells threshold needed to pass cos sim threshold below
params.condtn = "Symptomatic"

// - Comparison variable in metadata
// params.comps = "Symptomatic"

// - Knots number - should be estimated by tradeseq but generally 3-6 should be fine. So just using 4 here
params.knots = 4

// - number of genes threshold needed to be fitted GAM for
params.gns_cut = 5

// - tradeseq object output name
params.ts_out = "batch_ts.RDS"

// - condition test DE object output name
params.cond_de_out = "batch_de.RDS"

// - compute resources for first process
ncores="10"
mem="100 GB"

// - Standard scripts for running souporcell
params.scrpt = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/tradeseq_batch_nf.R"

// - Create channels 

// input_ch = Channel
// 				.fromPath(params.sce_obj)
//                 .collect()
// 				.map { file -> tuple(file[0].getParent().toString().split('\\/')[14], file[0], file[1], file[2], file[3], file[0].getParent(), file[0].getBaseName()) }

sce_obj_ch = Channel.fromPath("${params.w_dir}_sce.RDS").map { file -> tuple(file.getParent().toString().split('\\/')[14], file.getBaseName().toString().replace("_sce", ""), file ) }
ptime_obj_ch = Channel.fromPath("${params.w_dir}_ptime.RDS").map { file -> tuple(file.getParent().toString().split('\\/')[14], file.getBaseName().toString().replace("_ptime", ""), file ) }
cw_obj_ch = Channel.fromPath("${params.w_dir}_cw.RDS").map { file -> tuple(file.getParent().toString().split('\\/')[14], file.getBaseName().toString().replace("_cw", ""), file ) }
u_don_ch = Channel.fromPath("${params.w_dir}_u_don.RDS").map { file -> tuple(file.getParent().toString().split('\\/')[14], file.getBaseName().toString().replace("_u_don", ""), file ) }

input_ch = sce_obj_ch
            .combine(ptime_obj_ch, by: [0,1])
            .combine(cw_obj_ch, by: [0,1])
            .combine(u_don_ch, by: [0,1])
            .map { file -> tuple(file[0], file[1], file[2], file[3], file[4], file[5], file[5].getParent() ) }

input_ch.view()


// - Run souporcell assuming K=1 - remapping, variant calling and vartrix done only once and reused for other Ks
process TRADESEQ_DE {
    memory "${mem}"
    cpus "${ncores}"

    // errorStrategy 'ignore' // MSC66 giving error

    tag "tradeseq on ${sample_nm}"

    // publishDir params.outdir_index_temp, mode: "copy"
    publishDir "${out_p}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), val(obj_nm), path(sce_obj), path(ptime_obj), path(cw_obj), path(u_don_obj), val(out_p)
	
    output:
    tuple val(sample_nm), path("${obj_nm}_${params.ts_out}"), path("${obj_nm}_${params.cond_de_out}")

	script:
	"""
    module load HGI/softpack/groups/team222/Pf_scRNAseq_DE/5

    Rscript ${params.scrpt} ${sce_obj}  ${ptime_obj}  ${cw_obj} ${u_don_obj} ${params.condtn} ${params.knots} ${params.gns_cut} ${obj_nm}_${params.ts_out} ${obj_nm}_${params.cond_de_out}

	"""			
    
}


workflow {
    
    tseq_out_ch = TRADESEQ_DE(input_ch)
    tseq_out_ch.view()

}

