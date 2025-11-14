#!/usr/bin/env nextflow

// - sample id decode file
params.id_decode = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/irods_id_sk21_mdata_cln.csv" 

// - 2022 cellranger output
params.bam = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf/{5736STDY13782211,5736STDY13782212,5736STDY13782213,5736STDY13782214,5736STDY13782215, 5736STDY13437393}/outs/possorted_genome_bam.bam" 

// - cell barcodes
params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf/5736STDY*/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

// - output directory
params.o_dir= "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf"

// - compute resources for first process
ncores="20"
mem="120 GB"

// - souporcell output folder
params.soup_dir = "soupc_k5_verif"

// - Standard scripts for running souporcell
params.scrpt = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/soupc_v25_mnmap_hsat_310.sh"

// - Reference fasta
// params.ref  = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/hisat_refs/Pfalciparum.genome.fasta"

// - Create channels 
id_ch = Channel
				.fromPath(params.id_decode)
				.splitCsv( header: true )
				.map { file -> tuple( file.irods_id, file.sample_nm ) }
				// .map { file -> tuple( "${file.irods_id}, ${file.sample_nm}" ) } NB - this does not work

// id_ch.view()

bam_ch = Channel
				.fromPath(params.bam)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file, file+'.bai') }
				// .combine(id_ch, by:0)
				// .map { file -> tuple(file[3], file[1], file[2]) }
// bam_ch.view()

bcodes_ch = Channel
				.fromPath(params.bcodes)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file) }

// bcodes_ch.view()

bcodes_bam_ch = bcodes_ch
                    .combine(bam_ch, by:0)
                    .combine(id_ch, by:0)
                    .map { file -> tuple(file[4], file[1], file[2], file[3]) }


// - Algnment mapper tuples including references
hsat_tup = ['hsat', 'HISAT2', '/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/Pf66_hisat_refs/genome_w_tran_ref/genome_tran.fasta'] // reference generated using /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/scripts/hisat2_ref_build.sh
minmap_tup = ['minmap', 'minimap2', '/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/PlasmoDB-66_Pfalciparum3D7_Genome.fasta']

algn_ch = Channel.from([hsat_tup, minmap_tup])

bcodes_bam_al_ch = bcodes_bam_ch
                        .combine(algn_ch)

bcodes_bam_al_ch.view()


k_ch = Channel.from( 1..20 )
// k_ch.view()

// - Run souporcell assuming K=1 - remapping, variant calling and vartrix done only once and reused for other Ks
process SOUPC1 {
    memory "${mem}"
    cpus "${ncores}"

    // errorStrategy 'ignore'

    tag "Souporcell k1 ${sample_nm} ${algn_nm} reads"

    // publishDir params.outdir_index_temp, mode: "copy"
    publishDir "$params.o_dir/${sample_nm}/${params.soup_dir}/${algn_nm}/parent/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), path(bcodes), path(bam), path(bai), val(algn_nm), val(algnr), val(ref)
	
    output:
    tuple val(sample_nm), path(bcodes), path(bam), path(bai), val(algn_nm), val(algnr), val(ref), path(out_dir)

	script:
	"""
    
    ${params.scrpt} ${bam} ${bcodes} 5 ${algnr} ${ref} out_dir ${ncores}
    
	"""			
    
}

workflow {
    
    soupc1_ch = SOUPC1(bcodes_bam_al_ch)
    soupc1_ch.view()


}

