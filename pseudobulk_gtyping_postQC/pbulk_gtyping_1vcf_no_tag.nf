#!/usr/bin/env nextflow

/*
 * Pseudobulk genotyping of strain clusters from scRNAseq samples
 * - Performs variant calling and allele counting for clusters across multiple samples
 * - Uses FreeBayes, bcftools, and Vartrix
 * - Downsampling and dry-run options for testing
 */

// =======================
// PARAMETERS
// =======================

// Output directory for results
params.odir = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/all/"

// Directory for saving intermediate results (change as needed)
params.sv_dir = "pbulk_gtypes_cln_1VCF"

// Input BAM files (adjust glob as needed)
params.bam = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/${params.sv_dir}/strain_k1_*/new_bams/*_sset_sorted_rg.bam"

// Input cell barcode files
params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/${params.sv_dir}/bcodes/minmap/strain*/*.tsv"

// Reference FASTA
params.ref = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/PlasmoDB-66_Pfalciparum3D7_Genome.fasta"

// Downsampling BAMs for testing
params.downsample_bam = true
params.downsample_fraction = 01  // Fraction to downsample to (e.g. 01 = 1%)

// FreeBayes ploidy
params.ploidy = 1

// FreeBayes and bcftools parameters
fb_spec = Channel.value("-iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --min-alternate-fraction 0.2")
bcf_spec = Channel.value("view --max-alleles 2 -i 'QUAL>=75 && INFO/DP>=40'")

// Vartrix dry-run options
params.vartrix_dryrun = true
params.vartrix_dryrun_n = 10    // Number of barcodes to use in dry run

// =======================
// INPUT CHANNELS
// =======================

// Channel for BAM files, extracting metadata from file path
bam_ch = Channel
    .fromPath(params.bam)
    .map { file ->
        tuple(
            file.getParent().toString().split('/')[13], // sample name
            file.getParent().toString().split('/')[15].split('_')[-1], // alignment name - minimap or hisat2
            file.getParent().toString().split('/')[15].split('_minmap')[0].split('_k')[0], // group - strain or strain+stage
            file.baseName.toString().split('_sset_sorted_rg_don')[0], // strain stage combination eg MSC1_SC1_Asexual
            file, // bam file path
            file + '.bai' // bam file path index
        )
    }

//// bam_ch.view()

// Channel for barcode files, extracting metadata from file path
bcodes_ch = Channel
    .fromPath(params.bcodes)
    .map { file ->
        tuple(
            file.getParent().toString().split('/')[13], // sample name
            file.getParent().toString().split('/')[16], // alignment name
            file.getParent().name.split('_k')[0],       // group
            file.baseName,                              // strain stage
            file
        )
    }

//// bcodes_ch.view()

// bcodes_tagd_bam_ch = bcodes_ch
//                 .combine(bam_ch, by: [0,1,2,3])


// Reference FASTA and index
ref_ch = Channel.fromPath(params.ref)
    .map { file -> tuple(file, file + '.fai') }

// =======================
// PROCESS DEFINITIONS
// =======================

/*
 * SPLIT_GENOME
 * Split reference genome into regions for parallel variant calling
 */
process SPLIT_GENOME {
    memory '30 GB'
    cpus '1'

    tag "Split fasta index reference into 1 MB regions"

    input:
    tuple path(ref), path(ref_fai)

    output:
    path "regions.bed"

    script:
    """
   
    # Create 1 Mb genomic chunks (adjust size as needed)
    /software/team222/jr35/freebayes/fasta_generate_regions.py ${ref_fai} 1000000 > regions.bed
    """
}

/*
 * DOWNSAMPLE_BAM
 * Randomly downsample BAM files for testing
 */
process DOWNSAMPLE_BAM {
    memory '20 GB'
    cpus '1'
    queue 'normal'

    tag "Downsample the bam file for ${grp} pseudobulk cluster ${strn_stg} from ${sample_nm} aligned using ${algn_nm}"

    input:
    tuple val(sample_nm), val(algn_nm), val(grp), val(strn_stg), path(bam), path(bai), val(fraction)

    output:
    tuple val(sample_nm), val(algn_nm), val(grp), val(strn_stg), path("*.bam"), path("*.bai")

    script:
    """
    module load samtools/1.20--h50ea8bc_0
    samtools view -b -s 42.${fraction} ${bam} > \$(basename ${bam})_dsampled.bam
    samtools index \$(basename ${bam})_dsampled.bam
    """
}

/*
 * FREEBAYES
 * Variant calling for each region using FreeBayes and filtering with bcftools
 */
process FREEBAYES {
    memory '50 GB'
    cpus '2'
    queue 'normal'

    errorStrategy 'ignore'

    tag "Freebayes variant calling between ${grp} pseudobulk clusters for region ${region} for ${algn_nm}"

    input:
    tuple val(algn_nm), val(grp), path(bam_files), path(bai_files), path(ref), path(ref_fai), val(region)
    val fb_s
    val bcf_s

    output:
    tuple val(algn_nm), val(grp), path('*_biallelic_*.vcf.gz'), path('*_biallelic_*.vcf.gz.tbi'), path('*_multiallelic_*.vcf.gz'), path('*_multiallelic_*.vcf.gz.tbi')

    script:
    """
    module load samtools/1.20--h50ea8bc_0
    module load bcftools/1.20--h8b25389_0

    /software/team222/jr35/freebayes/freebayes_1.3.10 -f ${ref} ${fb_s} --ploidy ${params.ploidy} -r ${region} ${bam_files} | bgzip -c > fb_multiallelic_${region}.vcf.gz
    tabix -p vcf fb_multiallelic_${region}.vcf.gz

    bcftools ${bcf_s} fb_multiallelic_${region}.vcf.gz -Oz -o fb_biallelic_${region}.vcf.gz
    tabix -p vcf fb_biallelic_${region}.vcf.gz
    """
}

/*
 * MERGE_VCFS
 * Merge VCFs from all regions into a single file per sample/strain
 */
process MERGE_VCFS {
    memory '50 GB'
    cpus '1'
    queue 'normal'

    tag "Merge VCFs across regions for ${grp} pseudobulk clusters"

    publishDir "$params.odir/${params.sv_dir}/${grp}_kfinal_${algn_nm}/p${params.ploidy}/", mode: 'copy', pattern: "*.vcf", overwrite: true

    input:
    tuple val(algn_nm), val(grp), path(biallelic_vcfs), path(biallelic_vcfs_tbi), path(multiallelic_vcfs), path(multiallelic_vcfs_tbi)

    output:
    tuple val(algn_nm), val(grp), path("fb_biallelic.vcf"), path("fb_multiallelic.vcf")

    script:
    """
    module load bcftools/1.20--h8b25389_0

    bcftools concat -a -O v -o fb_biallelic.vcf ${biallelic_vcfs.join(' ')}
    # tabix -p vcf fb_multiallelic.vcf.gz

    bcftools concat -a -O v -o fb_multiallelic.vcf ${multiallelic_vcfs.join(' ')}
    # tabix -p vcf fb_multiallelic.vcf.gz
    """
}

/*
 * VARTRIX_UMI
 * Count alleles per barcode using Vartrix
 */
process VARTRIX_UMI {
    memory '100 GB'
    cpus '3'
    
    errorStrategy 'ignore'

    tag "Vartrix allele counting on ${grp} ${strn_stg} from ${sample_nm}"

    input:
    tuple val(sample_nm), val(algn_nm), val(grp), val(strn_stg), path(bcodes), path(bam), path(bai), path(bi_vcf)

    output:
    tuple val(sample_nm), val(algn_nm), val(grp), val(strn_stg), path ("${strn_stg}_*")

    script:
    """
    export PATH="\$PATH:/software/team222/jr35/vartrix/vartrix-1.1.22/target/release"

    # Use subset of barcodes for dry run if enabled
    if [ "${params.vartrix_dryrun}" = "true" ]; then
        echo "Running Vartrix dry run: using first ${params.vartrix_dryrun_n} barcodes"
        head -n ${params.vartrix_dryrun_n} ${bcodes} > subset_bcodes.tsv
        use_bcodes=subset_bcodes.tsv
    else
        use_bcodes=${bcodes}
    fi

    vartrix --umi --out-variants ${strn_stg}_variants.txt --mapq 30 -b ${bam} -c \${use_bcodes} --scoring-method coverage --ref-matrix ${strn_stg}_ref.mtx --out-matrix ${strn_stg}_alt.mtx -v ${bi_vcf} --fasta ${params.ref}
    """
}

/*
 * VARTRIX_COPY
 * Organize Vartrix output files for downstream analysis
 */
process VARTRIX_COPY {
    memory '40 GB'
    cpus '1'
    debug true

    tag "Copy vartrix output for ${grp}"

    publishDir "$params.odir/${params.sv_dir}/${grp}_kfinal_${algn_nm}/p${params.ploidy}/", mode: 'copy', pattern: "vartrix_biallelic", overwrite: true

    input:
    tuple val(algn_nm), val(grp), path(vtx_files)

    output:
    tuple val(algn_nm), val(grp), path ("vartrix_biallelic")

    script:
    """
    mkdir vartrix_biallelic
    mv ${vtx_files} vartrix_biallelic/
    """
}

// =======================
// WORKFLOW DEFINITION
// =======================

workflow {

    // Split genome into regions for parallel variant calling
    SPLIT_GENOME(ref_ch)
    regions_ch = SPLIT_GENOME.out
        .splitText()
        .map { it.trim() }
        .filter { it }

    //regions_ch.view()

    // Downsample BAMs if enabled
    def downsample = params.downsample_bam.toString().toLowerCase() in ['true', '1', 'yes']
    downsampled_ch = downsample ?
        bam_ch.map { tuple(it[0], it[1], it[2], it[3], it[4], it[5], params.downsample_fraction) } | DOWNSAMPLE_BAM
        :
        bam_ch

    // Group BAMs for each alignment/group
    tagbam_ch_all = downsampled_ch
        .map { file -> tuple(file[1], file[2], file[4], file[5]) }
        .groupTuple(by: [0,1])

    // Combine BAMs, reference, and regions for variant calling
    bam_region_ch = tagbam_ch_all
        .combine(ref_ch)
        .combine(regions_ch)

    // Run FreeBayes and filter variants to retain biallelic sites for vartrix and downstream analysis
    fb_vcf_ch = FREEBAYES(bam_region_ch, fb_spec, bcf_spec)
    fb_vcf_grouped_ch = fb_vcf_ch.groupTuple(by: [0, 1])

    // Merge VCFs from all regions
    fb_vcf_merge_ch = MERGE_VCFS(fb_vcf_grouped_ch)

    // Prepare VCFs for Vartrix
    fb_vcf_bi_ch = fb_vcf_merge_ch.map { file -> tuple(file[2], file[0], file[1]) }

    // Match barcodes to BAMs and VCFs
    bcodes_tagd_bam_ch = bcodes_ch.combine(downsampled_ch, by: [0,1,2,3])
    bcodes_fb_vcf_ch = bcodes_tagd_bam_ch
        .combine(fb_vcf_bi_ch, by: [1,2])
        .map { file -> tuple(file[2], file[0], file[1], file[3], file[4], file[5], file[6], file[7]) }

    // Run Vartrix allele counting
    vartx_ch = VARTRIX_UMI(bcodes_fb_vcf_ch)

    // Organize and group Vartrix outputs
    vartx_clean_ch = vartx_ch
        .map { file -> tuple(file[0], file[1], file[2], file[4]) }
        .groupTuple(by: [0, 1, 2])
        .map { gt ->
            def keys = gt.take(3)
            def values = gt.drop(3).flatten().collect()
            return keys + [values]
        }
        .map { file -> tuple(file[1], file[2], file[3]) }
        .groupTuple(by: [0,1])
        .map { file -> tuple(file[0], file[1], file[2].flatten().collect()) }

    // Publish Vartrix outputs 
    vartx_cpy_ch = VARTRIX_COPY(vartx_clean_ch)
}
