#!/usr/bin/env nextflow

// - Directory to store output - need to differentiate between preQC genotyping so as to verify souporcell assignments and postQC genotyping so as to perform relatedness analysis etc
// params.sv_dir = "pbulk_gtypes_preQC"
// params.sv_dir = "pbulk_gtypes_GE_postQC"
params.sv_dir = "pbulk_gtypes_postQC_cln"

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
params.ref  = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/PlasmoDB-66_Pfalciparum3D7_Genome.fasta"

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
// bcodes_bam_ch.view()

// - Processes
// - Subset minimap bam file from souporcell based on the barcodes in each pseudobulk group using subset bam from 10X

process SUBSET {
    memory '60 GB'
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

    tag "Tagging ${strn_stg}  reads"

    // publishDir "$params.odir/${sample_nm}/${params.sv_dir}/${grp}_${algn_nm}/nw_bams", mode: 'copy', pattern: "*sset_sorted_rg.bam*"
        
    input:
    tuple val(sample_nm), val(algn_nm), val(grp), val(strn_stg), path(bcodes), path(bam), path(bai)

    output:
    tuple val(sample_nm), val(algn_nm), val(grp), val(strn_stg), path(bcodes), path('*_sset_sorted_rg.bam'), path('*_sset_sorted_rg.bam.bai')
    
    script:
    """
    module load samtools/1.20--h50ea8bc_0

    samtools addreplacerg -r ID:${strn_stg} -r SM:${strn_stg} $bam -o ${strn_stg}_sset_sorted_rg.bam
    samtools index ${strn_stg}_sset_sorted_rg.bam
    """
}

// variant calling variables
// params.ref  = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/PlasmoDB-52_Pfalciparum3D7_Genome.fasta"

ref_ch = Channel.value(params.ref)
params.ploidy=(1)

// freebayes parameters
// fb_spec=Channel.of("-iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --min-alternate-fraction 0.2" "-C 2 -q 20 -n 3 -E 3 -m 30 --min-coverage 6 --min-alternate-fraction 0.2 --theta 0.01")
fb_spec=Channel.value("-iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --min-alternate-fraction 0.2")

// Bcftools vcf processing parameters
// bcf_spec=Channel.of("view --max-alleles 2" "norm --multiallelics -")
bcf_spec=Channel.value("view --max-alleles 2")

// Run FreeBayes with the specified parameters using ploidy of 1 and other parameters specified above

process FREEBAYES {
    memory '100 GB'
    cpus '8'

    errorStrategy 'ignore'

    tag "Freebayes variant calling on ${grp} ${sample_nm}"

    publishDir "$params.odir/${sample_nm}/${params.sv_dir}/${grp}_${algn_nm}/p${params.ploidy}/", mode: 'copy', pattern: "*.vcf", overwrite: true
    
    input:
    tuple val(sample_nm), val(algn_nm), val(grp), path(bam), path(bai)
    val fb_s
    val bcf_s

    output:
    tuple val(sample_nm), val(algn_nm), val(grp), path('*_biallelic.vcf'), path('*_multiallelic.vcf')

    script:
    """
    module load samtools/1.20--h50ea8bc_0
    module load HGI/softpack/groups/team222/SNP_analysis/1
    module load bcftools/1.20--h8b25389_0

    freebayes -f $params.ref $fb_s --ploidy $params.ploidy $bam > fb_multiallelic.vcf 
    bcftools $bcf_s fb_multiallelic.vcf -o fb_biallelic.vcf
    """
}

process VARTRIX_UMI {
    memory '100 GB'
    cpus '8'
    
    errorStrategy 'ignore'

    tag "Vartrix allele counting on ${strn_stg} from ${grp} ${sample_nm}"
        
    input:
    tuple val(sample_nm), val(algn_nm), val(grp), val(strn_stg), path(bcodes), path(bam), path(bai), path(bi_vcf)
    
    output:
    tuple val(sample_nm), val(algn_nm), val(grp), val(strn_stg), path ("${strn_stg}_*")
    
    script:
    """
    export PATH="\$PATH:/software/team222/jr35/vartrix/vartrix-1.1.22/target/release"

    vartrix --umi --out-variants ${strn_stg}_variants.txt --mapq 30 -b ${bam} -c ${bcodes} --scoring-method coverage --ref-matrix ${strn_stg}_ref.mtx --out-matrix ${strn_stg}_alt.mtx -v ${bi_vcf} --fasta ${params.ref}
    
    """
}

process VARTRIX_COPY {
    memory '40 GB'
    cpus '1'
    debug true

    tag "Vartrix copy from ${grp} ${sample_nm}"
    
    publishDir "$params.odir/${sample_nm}/${params.sv_dir}/${grp}_${algn_nm}/p${params.ploidy}/", mode: 'copy', pattern: "vartrix_biallelic", overwrite: true
    
    input:
    tuple val(sample_nm), val(algn_nm), val(grp), path(vtx_files)
    
    output:
    tuple val(sample_nm), val(algn_nm), val(grp), path(vartrix_biallelic)
    
    script:
    """
    mkdir vartrix_biallelic
    mv ${vtx_files} vartrix_biallelic/
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

    fb_vcf_ch = FREEBAYES(tagbam_ch_all,  fb_spec, bcf_spec)

    // fb_vcf_ch.view()
    fb_vcf_bi_ch = fb_vcf_ch.map { file -> tuple(file[0], file[1], file[2], file[3]) }
    // fb_vcf_bi_ch.view()

    bcodes_fb_vcf_ch = tagd_bams_ch
                        .combine(fb_vcf_bi_ch, by: [0, 1, 2])

    vartx_ch = VARTRIX_UMI(bcodes_fb_vcf_ch)
                    .map { file -> tuple(file[0], file[1], file[2], file[4]) }
                    .groupTuple(by: [0, 1, 2])
                    .map { gt ->
                            def keys = gt.take(3)
                            def values = gt.drop(3).flatten().collect()
                            return keys + [values]
                        }

    vartx_ch.view()

    vartx_cpy_ch =VARTRIX_COPY(vartx_ch)

    vartx_cpy_ch.view()
            
}
