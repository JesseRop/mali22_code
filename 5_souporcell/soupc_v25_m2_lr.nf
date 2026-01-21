#!/usr/bin/env nextflow

// NOTE - This script only runs the entire souporcell pipeline for the first K and then uses the bam, vcf and vartrix filed generated from this for other Ks to save on resources. 
// The shortcut gives identical results to those from running the pipeline from scratch as verified by output from the K=5 rerun using soupc_v25_m2_k5_verif.nf script

// PART 1 - VARIABLES AND CHANNELS

// // - 2022 long read bam with correct cell barcodes
params.bam = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf/5736STDY134373{63,72,89}/outs/possorted_genome_bam.bam" 


// - filtered cell barcodes - for first pass souporcell
params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/Pf/MSC*/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

// - output directory to store souporcell output
params.o_dir= "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/"

// souporcell output folder - this inside the directory above and we may run different iterations of souporcell hence each should go to a different directory
params.soup_dir = "soupc"


// - compute resources for first process
ncores="10"
mem="100 GB"


// - Standard scripts for running souporcell - this is put in a different scripts since it is used for other projects
params.scrpt = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/soupc_v25_mnmap_hsat_310_lr.sh"

// - Reference fasta
params.ref_file  = "PlasmoDB-66_Pfalciparum3D7_Genome.fasta"

// tag for UMI - long read and short read  uses a different tag
params.tag  = "XM"

// Bam files channel for different donors
bam_ch = Channel
				.fromPath(params.bam)
				// .map { file -> tuple(file.getParent().toString().split('\\/')[13], file, file+'.bai') }
                .map { file -> tuple(file.baseName.split('_')[0], file, file+'.bai') }
				// .combine(id_ch, by:0)
				// .map { file -> tuple(file[3], file[1], file[2]) }
bam_ch.view()

// Barcode channel for different donors
bcodes_ch = Channel
				.fromPath(params.bcodes)
				.map { file -> tuple(file.baseName.split('_')[0], file) }

bcodes_ch.view()

// Combine barcodes and bam file for each donor
bcodes_bam_ch = bcodes_ch
                    .combine(bam_ch, by:0)

bcodes_bam_ch.view()

// // - Algnment mapper tuples including name for folder, name needed by souporcell and references path
minmap_tup = ["minmap", "minimap2", "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/${params.ref_file}"]

// Channel for the alignment
algn_ch = Channel.from([minmap_tup])

// Combine barcodes and bam file channel with alignment channel for each donor
bcodes_bam_al_ch = bcodes_bam_ch
                        .combine(algn_ch)

bcodes_bam_al_ch.view()

// Channel for the souporcell range of strain numbers
k_ch = Channel.from( 1..20 )
// // k_ch.view()

// PART 2 - PROCESSES

// NOTE - to ensure the variables and channels above are working as expected comment out everything below this line, save and then run and inspect the output. 
// Once happy with the paths and variables in the channels then undo the commenting and proceed to run


// - Run souporcell assuming K=1 - remapping, variant calling and vartrix done only once and reused for other Ks
process SOUPC1 {
    memory "${mem}"
    cpus "${ncores}"

    // errorStrategy 'ignore' // MSC66 giving error

    tag "Souporcell k1 ${sample_nm} ${algn_nm} reads"

    publishDir "$params.o_dir/${sample_nm}/${params.soup_dir}/${algn_nm}/parent/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_nm), path(bcodes), path(bam), path(bai), val(algn_nm), val(algnr), val(ref)
	
    output:
    tuple val(sample_nm), path(bcodes), path(bam), path(bai), val(algn_nm), val(algnr), val(ref), path(out_dir)

	script:
	"""
    
    ${params.scrpt} ${bam} ${bcodes} 1 ${algnr} ${ref} out_dir ${ncores} ${params.tag}
    
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
	${params.scrpt} ${bam} ${bcodes} ${k} ${algnr} ${ref} out_dir_k${k} ${ncores}  ${params.tag}

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
    tuple val(sample_nm), path('clusters_log_likelihoods.txt')
    
	script:
	"""
    grep "" $clst_err > clusters_log_likelihoods.txt
    
	"""			
}

// PART 2 - WORKFLOWS

workflow {
    
    soupc1_ch = SOUPC1(bcodes_bam_al_ch)
    // soupc1_ch.view()
    
    soupc2plus_ch = SOUPC_2PLUS(soupc1_ch, k_ch)
    // soupc2plus_ch.view()
    
    soupc2plus_ll_ch = soupc2plus_ch
        .map { file -> tuple(file[0], file[1], file[5]) }
        .groupTuple(by: [0,1])
    
    soupc2plus_ll_ch.view()

    LL_KNEE_PLOT(soupc2plus_ll_ch)

}

