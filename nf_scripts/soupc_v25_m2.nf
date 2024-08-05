#!/usr/bin/env nextflow

// NOTE - This script only runs the entire souporcell pipeline for the first K and then uses the bam, vcf and vartrix filed generated from this for other Ks to save on resources. 
// The shortcut gives identical results to those from running the pipeline from scratch as verified by output from the K=5 rerun using soupc_v25_m2_k5_verif.nf script

// - sample id decode file. second file includes mixed infections
// params.id_decode = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/irods_id_sk21_mdata_cln.csv"
params.id_decode = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/pf_solo_mixed_mdata_decode.csv"  

// - 2022 cellranger output - for all infections
// params.bam = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf/5736STDY*/outs/possorted_genome_bam.bam" 

// // - 2022 cellranger output - for mixed infections
params.bam = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf/5736STDY134373{63,72,89}/outs/possorted_genome_bam.bam" 


// - Cellranger filtered cell barcodes - for first pass souporcell
params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/Pf/MSC*/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

// - NB!!! -DO NOT DELETE
// - Clean cell barcodes - use this for the second pass - limiting souporcell only to the QCd (outlier low/high count and doublets removed) barcodes
// params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/soupc_GE_postQC/barcodes.tsv.gz"

// - output directory
params.o_dir= "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/"

// - compute resources for first process
ncores="15"
mem="120 GB"

// - Cellranger filtered cell barcodes for first pass souporcell
// souporcell output folder
params.soup_dir = "soupc"

// - Clean cell barcodes - use this for the second pass
// souporcell output folder
// params.soup_dir = "soupc_GE_postQC"

// - Standard scripts for running souporcell
params.scrpt = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/soupc_v25_mnmap_hsat_310.sh"

// - Reference fasta
params.ref_file  = "PlasmoDB-66_Pfalciparum3D7_Genome.fasta"
params.hsat_ref_file  = "Pf66"

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

// Old channel script for when the barcodes file was being obtained from sk21's directory 
// bcodes_bam_ch = bcodes_ch
//                     .combine(bam_ch, by:0)
//                     .combine(id_ch, by:0)
//                     .map { file -> tuple(file[4], file[1], file[2], file[3]) }

bcodes_bam_ch = id_ch
                    .combine(bam_ch, by:0)
                    .map { file -> tuple(file[1], file[0], file[2], file[3]) }
                    .combine(bcodes_ch, by:0)
                    .map { file -> tuple(file[0], file[4], file[2], file[3]) }

// bcodes_bam_ch.view()

// - Algnment mapper tuples including references
hsat_tup = ["hsat", "HISAT2", "/lustre/scratch126/tol/teams/lawniczak/users/jr35/genomes_gtfs_indexs/${params.hsat_ref_file}_hisat_refs/genome_w_tran_ref/genome_tran.fasta"] // reference generated using /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/scripts/hisat2_ref_build.sh
minmap_tup = ["minmap", "minimap2", "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/${params.ref_file}"]
// lr_tup = ["lr_ref", "minimap2", "/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/sk21_long_read_ref/Ref_pxPlaFalc47_hum_on_PfDB66_no_contigs_no_mit.fa"]


// algn_ch = Channel.from([hsat_tup, minmap_tup, lr_tup])
algn_ch = Channel.from([hsat_tup, minmap_tup])
// algn_ch = Channel.from([minmap_tup])

bcodes_bam_al_ch = bcodes_bam_ch
                        .combine(algn_ch)

bcodes_bam_al_ch.view()


k_ch = Channel.from( 1..20 )
// k_ch.view()

// - Run souporcell assuming K=1 - remapping, variant calling and vartrix done only once and reused for other Ks
process SOUPC1 {
    memory "${mem}"
    cpus "${ncores}"

    // errorStrategy 'ignore' // MSC66 giving error

    tag "Souporcell k1 ${sample_nm} ${algn_nm} reads"

    // publishDir params.outdir_index_temp, mode: "copy"
    publishDir "$params.o_dir/${sample_nm}/${params.soup_dir}/${algn_nm}/parent/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), path(bcodes), path(bam), path(bai), val(algn_nm), val(algnr), val(ref)
	
    output:
    tuple val(sample_nm), path(bcodes), path(bam), path(bai), val(algn_nm), val(algnr), val(ref), path(out_dir)

	script:
	"""
    
    ${params.scrpt} ${bam} ${bcodes} 1 ${algnr} ${ref} out_dir ${ncores}
    
	"""			
    
}


// - Run souporcell assuming K=2..etc - reuse bam, vcf and mtx generated for K=1
process SOUPC_2PLUS {
    memory '100 GB'
    cpus '5'

    errorStrategy 'ignore'

    tag "Souporcell k${k} ${sample_nm} ${algn_nm} reads"

    publishDir "$params.o_dir/${sample_nm}/${params.soup_dir}/${algn_nm}/k${k}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), path(bcodes), path(bam), path(bai), val(algn_nm), val(algnr), val(ref), path(out_dir)
    each(k)
		
    output:
    tuple val(sample_nm), val(algn_nm), path('clusters.tsv'), path('cluster_genotypes.vcf'), path('ambient_rna.txt'), path('clusters_log_likelihoods*')  

	script:
	"""
    
    mkdir out_dir_k${k}
    cp ${out_dir}/{vartrix.done,alt.mtx,ref.mtx,barcodes.tsv,bcftools.err,souporcell_merged_sorted_vcf.vcf.gz,souporcell_merged_sorted_vcf.vcf.gz.tbi,retagging.done,souporcell_minimap_tagged_sorted.bam.bai,souporcell_minimap_tagged_sorted.bam,retag.err,remapping.done,minimap.err,fastqs.done} out_dir_k${k}/
    echo 'out_dir_k${k}/souporcell_merged_sorted_vcf.vcf.gz' > out_dir_k${k}/variants.done
	${params.scrpt} ${bam} ${bcodes} ${k} ${algnr} ${ref} out_dir_k${k} ${ncores}

    cp out_dir_k*/cluster* out_dir_k*/ambient_rna.txt .

    grep -H 'best total log probability' out_dir_k${k}/clusters.err > clusters_log_likelihoods_${k}.txt

    rm out_dir_k${k}/{vartrix.done,alt.mtx,ref.mtx,barcodes.tsv,bcftools.err,souporcell_merged_sorted_vcf.vcf.gz,souporcell_merged_sorted_vcf.vcf.gz.tbi,retagging.done,souporcell_minimap_tagged_sorted.bam.bai,souporcell_minimap_tagged_sorted.bam,retag.err,remapping.done,minimap.err,fastqs.done}
    
	"""			
}

// - Put together the log likelihood for knee plot
process LL_KNEE_PLOT {
    memory '10 GB'
    cpus '1'

    errorStrategy 'ignore'

    tag "Souporcell log likelihood values ${sample_nm} reads"

    publishDir "$params.o_dir/${sample_nm}/${params.soup_dir}/${algn_nm}/", mode: 'copy', overwrite: true

    input:
	tuple val(sample_nm), val(algn_nm), path(clst_err)
		
    output:
    // tuple val(sample_nm), path('*')
    tuple val(sample_nm), path('clusters_log_likelihoods.txt')
    
	script:
	"""
    grep "" $clst_err > clusters_log_likelihoods.txt
    
	"""			
}

workflow {
    
    soupc1_ch = SOUPC1(bcodes_bam_al_ch)
    // soupc1_ch.view()

    // bcodes_bam_soupc_ch = bcodes_bam_al_ch
    //     .map { file -> tuple(file[0], file[1], file[2], file[3], file[5], file[6]) }
    //     .combine(soupc1_ch, by:[0,1])
    // bcodes_bam_soupc_ch.view()
    
    soupc2plus_ch = SOUPC_2PLUS(soupc1_ch, k_ch)
    // soupc2plus_ch = SOUPC_2PLUS(bcodes_bam_soupc_ch, k_ch)
    // soupc2plus_ch.view()
    soupc2plus_ll_ch = soupc2plus_ch
        .map { file -> tuple(file[0], file[1], file[5]) }
        .groupTuple(by: [0,1])
    
    soupc2plus_ll_ch.view()

    LL_KNEE_PLOT(soupc2plus_ll_ch)

}

