#!/usr/bin/env nextflow

// - Directory to store output - need to differentiate between preQC genotyping so as to verify souporcell assignments and postQC genotyping so as to perform relatedness analysis etc
// params.sv_dir = "pbulk_gtypes_preQC"
// params.sv_dir = "pbulk_gtypes_GE_postQC"
params.sv_dir = "pbulk_gtypes_postQC_cln"

// Whether to publish subsetted bam files
params.save_bam_sset = false

// - 2022 cellranger output
params.bam = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/soupc/*/parent/possorted_genome_bam.bam" 

// - cell barcodes
// params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC50_SLO/${params.sv_dir}/bcodes/*/stage_afm_strain_k[0-9]*/*.tsv"
params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/${params.sv_dir}/bcodes/minmap/*/*.tsv"

// - output directory
params.odir= "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/"

// ncores=10
// params.mapper="minmap"

// - Standard scripts for running souporcell
params.scrpt = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/cr_subset_bam_linux.sh"

// - Reference fasta
// params.ref  = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/PlasmoDB-66_Pfalciparum3D7_Genome.fasta"

// params.sv_dir = "pbulk_gtypes_postQC"

// - Create channels
bam_ch = Channel
				.fromPath(params.bam)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file.getParent().toString().split('\\/')[15], file, file+'.bai') }

// bam_ch.view()

bcodes_ch = Channel
				.fromPath(params.bcodes)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file.getParent().toString().split('\\/')[16], file.getParent().name, file.baseName, file) }

// bcodes_ch.view()

bcodes_bam_ch = bcodes_ch.combine(bam_ch, by: [0,1])
bcodes_bam_ch.view()

// - Processes
// - Subset minimap bam file from souporcell based on the barcodes in each pseudobulk group using subset bam from 10X

process SUBSET {
    memory '100 GB'
    cpus '5'

    tag "Subset ${strn_stg}  reads"
    
    input:
    tuple val(sample_nm), val(algn_nm), val(grp), val(strn_stg), path(bcodes), path(bam), path(bai)
	path scrpt
	
    output:
    tuple val(sample_nm), val(algn_nm), val(grp), val(strn_stg), path(bcodes), path('*_sset_sorted.bam'), path('*_sset_sorted.bam.bai')

	script:
	"""
    ${params.scrpt} ${strn_stg} ${bcodes} ./ ${bam}
    
	"""			
    
}

// Tag each subsetted bam file using the pseudobulk group ID in preparation for freebayes genotyping
process TAGBAM {
    memory '40 GB'
    cpus '3'

    errorStrategy 'ignore' // Some error coming up "terminated for an unknown reason -- Likely it has been terminated by the external system"

    tag "Tagging ${sample_nm}_${strn_stg}  reads"

    // publish subsetted bam file for IGV etc only for the final output
    publishDir "$params.odir/${sample_nm}/${params.sv_dir}/${grp}_${algn_nm}/new_bams", mode: 'copy', pattern: "*_sset_sorted_rg_don.bam*", overwrite: true, enabled: params.save_bam_sset

    input:
    tuple val(sample_nm), val(algn_nm), val(grp), val(strn_stg), path(bcodes), path(bam), path(bai)

    output:
    tuple val(sample_nm), val(algn_nm), val(grp), val(strn_stg), path(bcodes), path('*_sset_sorted_rg_don.bam'), path('*_sset_sorted_rg_don.bam.bai')
    
    script:
    """
    module load samtools/1.20--h50ea8bc_0

    samtools addreplacerg -r ID:${sample_nm}_${strn_stg} -r SM:${sample_nm}_${strn_stg} $bam -o ${sample_nm}_${strn_stg}_sset_sorted_rg_don.bam
    samtools index ${sample_nm}_${strn_stg}_sset_sorted_rg_don.bam

    """
}


// - Workflow
workflow {
    
    sset_bams_ch = SUBSET(bcodes_bam_ch, params.scrpt)
    // sset_bams_ch.view()

    tagd_bams_ch = TAGBAM(sset_bams_ch)
    // tagd_bams_ch.view()

    tagbam_ch_all = tagd_bams_ch
                        .map { file -> tuple(file[0], file[1], file[2], file[5], file[6]) }
                        .groupTuple(by: [0, 1, 2])
    
    tagbam_ch_all.view()

            
}
