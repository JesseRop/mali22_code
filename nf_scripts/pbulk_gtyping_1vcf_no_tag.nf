#!/usr/bin/env nextflow

// - Directory to store output - need to differentiate between preQC genotyping so as to verify souporcell assignments and postQC genotyping so as to perform relatedness analysis etc
// params.sv_dir = "pbulk_gtypes_preQC"
// params.sv_dir = "pbulk_gtypes_GE_postQC"
// params.sv_dir = "pbulk_gtypes_preQC_cln"
params.sv_dir = "pbulk_gtypes_cln_1VCF"

// Whether to publish subsetted bam files
// params.save_bam_sset = false

// - 2022 cellranger output
// params.bam = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/soupc/*/parent/possorted_genome_bam.bam"
params.bam = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/${params.sv_dir}/strain_k1_*/new_bams/*_sset_sorted_rg.bam"  
			  
// - cell barcodes
// params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC50_SLO/${params.sv_dir}/bcodes/*/stage_afm_strain_k[0-9]*/*.tsv"
params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/${params.sv_dir}/bcodes/minmap/strain*/*.tsv"

// - output directory
params.odir= "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/all/"

// ncores=10
// params.mapper="minmap"

// - Standard scripts for running souporcell
// params.scrpt = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/cr_subset_bam_linux.sh"

// - Reference fasta
params.ref  = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/PlasmoDB-66_Pfalciparum3D7_Genome.fasta"

// params.sv_dir = "pbulk_gtypes_postQC"

// - Create channels
bam_ch = Channel
				.fromPath(params.bam)
				// .map { file -> tuple(file.getParent().toString().split('\\/')[13], file.getParent().toString().split('\\/')[15].split('_')[-1], file.getParent().toString().split('\\/')[15].split('_minmap')[0].split('_k')[0], file.baseName.toString().split('_sset_sorted_rg_don')[0].replaceFirst(/^MSC\d{2}_/, ''), file, file+'.bai') }
                .map { file -> tuple(file.getParent().toString().split('\\/')[13], file.getParent().toString().split('\\/')[15].split('_')[-1], file.getParent().toString().split('\\/')[15].split('_minmap')[0].split('_k')[0], file.baseName.toString().split('_sset_sorted_rg_don')[0], file, file+'.bai') }

// bam_ch.view()
// fb_vcf_bi_ch = Channel
// 				.fromPath("/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/pbulk_gtypes_preQC_cln/minmap/p1/fb_biallelic.vcf")
//                 .map { file -> tuple(file, "minmap")}
                
// fb_vcf_bi_ch.view()

bcodes_ch = Channel
				.fromPath(params.bcodes)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file.getParent().toString().split('\\/')[16], file.getParent().name.split('_k')[0], file.baseName, file) }
                

// bcodes_ch.view()

bcodes_tagd_bam_ch = bcodes_ch
                .combine(bam_ch, by: [0,1,2,3])

// bcodes_tagd_bam_ch.view()

tagbam_ch_all = bcodes_tagd_bam_ch
                        .map { file -> tuple(file[1], file[2], file[5], file[6]) }
                        .groupTuple(by: [0,1]) 

						
tagbam_ch_all.view()


ref_ch = Channel.value(params.ref)
params.ploidy=(1)

// freebayes parameters
// fb_spec=Channel.of("-iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --min-alternate-fraction 0.2" "-C 2 -q 20 -n 3 -E 3 -m 30 --min-coverage 6 --min-alternate-fraction 0.2 --theta 0.01")
fb_spec=Channel.value("-iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --min-alternate-fraction 0.2")

// Bcftools vcf processing parameters
// bcf_spec=Channel.of("view --max-alleles 2" "norm --multiallelics -")
bcf_spec=Channel.value("view --max-alleles 2 -i 'QUAL>=75 && INFO/DP>=40'") // Retaining variants that have overall quality/depth across all samples of 75/40

// Run FreeBayes with the specified parameters using ploidy of 1 and other parameters specified above


// - Processes

process FREEBAYES {
    memory '700 GB'
    cpus '30'
    queue 'hugemem'
    // conda '/software/team222/jr35/miniconda3/envs/freebayes_para'


    // errorStrategy 'ignore'

    tag "Freebayes variant calling on ${params.sv_dir}"

    // publishDir "$params.odir/${sample_nm}/${params.sv_dir}/${grp}_${algn_nm}/p${params.ploidy}/", mode: 'copy', pattern: "*.vcf", overwrite: true
    publishDir "$params.odir/${params.sv_dir}/${strn_stg}_kfinal_${algn_nm}/p${params.ploidy}/", mode: 'copy', pattern: "*.vcf", overwrite: true

    input:
    tuple val(algn_nm), val(strn_stg), path(bam), path(bai)
    val fb_s
    val bcf_s

    output:
    tuple val(algn_nm), val(strn_stg), path('*_biallelic.vcf'), path('*_multiallelic.vcf')

    script:
    """
    module load samtools/1.20--h50ea8bc_0
    module load bcftools/1.20--h8b25389_0
    # module load HGI/softpack/groups/team222/SNP_analysis/2
    # module load ISG/singularity/3.11.4
    # module load HGI/softpack/groups/team222/Variant_calling/2
    module load HGI/common/freebayes/v1.2.0 


    
    freebayes -f $params.ref $fb_s --ploidy $params.ploidy $bam > fb_multiallelic.vcf
    #/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code/bash_scripts/freebayes-parallel /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code/python_nbooks/pf_region_68.txt 48 -f $params.ref $fb_s --ploidy $params.ploidy $bam > fb_multiallelic.vcf
    bcftools $bcf_s fb_multiallelic.vcf -o fb_biallelic.vcf
    """
}

process VARTRIX_UMI {
    memory '150 GB'
    cpus '10'
    
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

    tag "Vartrix copy"
    
    publishDir "$params.odir/${params.sv_dir}/${strn_stg}_kfinal_${algn_nm}/p${params.ploidy}/", mode: 'copy', pattern: "vartrix_biallelic", overwrite: true
                
    input:
    tuple val(algn_nm), val(strn_stg), path(vtx_files)
    
    output:
    tuple val(algn_nm), val(strn_stg), path ("vartrix_biallelic")
    
    script:
    """
    mkdir vartrix_biallelic
    mv ${vtx_files} vartrix_biallelic/
    """
}

// - Workflow
workflow {
    
    fb_vcf_ch = FREEBAYES(tagbam_ch_all,  fb_spec, bcf_spec)

    // fb_vcf_ch.view()
    fb_vcf_bi_ch = fb_vcf_ch.map { file -> tuple(file[2], file[0], file[1]) }
    fb_vcf_bi_ch.view()

    bcodes_fb_vcf_ch = bcodes_tagd_bam_ch
                        // .combine(fb_vcf_bi_ch, by: [0, 1, 2])
                        .combine(fb_vcf_bi_ch, by: [1,2])
                        .map { file -> tuple(file[2], file[0], file[1], file[3], file[4], file[5], file[6], file[7]) }
                        // .join(fb_vcf_bi_ch, by: [1])

	// bcodes_fb_vcf_ch.view()

    vartx_ch = VARTRIX_UMI(bcodes_fb_vcf_ch)
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
                    

    // vartx_ch.view()

    vartx_cpy_ch =VARTRIX_COPY(vartx_ch)

    vartx_cpy_ch.view()
            
}
