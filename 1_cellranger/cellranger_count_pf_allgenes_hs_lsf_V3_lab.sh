#!/bin/bash
#BSUB -q long
#BSUB -G team222
#BSUB -J crc_pf[1-4]
#BSUB -cwd /lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/v3_atlas/data/cellranger_runs/Pf_all_genes/
#BSUB -o logs/cellranger_count_PlasmoDB68_all_genes_GRCh38_Pf_JC_combined_%J_%I.out
#BSUB -e logs/cellranger_count_PlasmoDB68_all_genes_GRCh38_Pf_JC_combined_%J_%I.err
#BSUB -n 20
#BSUB -R "select[mem>84000] rusage[mem=84000]"
#BSUB -M84000

core=20
mem=80

cellranger=/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/resources/software/cellranger/cellranger-9.0.1/bin/cellranger

mcadir=/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/
fastqdir=${mcadir}mali_field_runs/raw_data/10X_short_read_fastqs/36012
# outdir=${mcadir}mali_field_runs/v3_atlas/data/cellranger_runs/Pf_all_genes
#cr_ref=${mcadir}resources/genome_transcript_refs/PlasmoDB68_GRCh38.p14_set/cellranger/Pf3D7_PDB68
cr_ref=${mcadir}resources/genome_transcript_refs/PlasmoDB68_GRCh38.p14_set/cellranger/Pf3D7_PDB68_all_genes_GRCh38

## NOTE!! - Using a TSV file here to avoid issues with commas within quoted fields in the CSV
tsv=/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/v3_atlas/data/raw/v3_atlas_decode_cellranger.tsv
sample_id=$(awk -F'\t' -v idx=$LSB_JOBINDEX 'NR==idx+1 {print $2}' "$tsv")
id=$(awk -F'\t' -v idx=$LSB_JOBINDEX 'NR==idx+1 {print $2}' "$tsv")

echo "Processing sample $LSB_JOBINDEX:"
echo "Sample ID: $sample_id"
echo "Output ID: $id"
echo "-------------------"


${cellranger} count --id=$id \
                 --transcriptome=$cr_ref \
                 --fastqs=$fastqdir \
                 --sample=$sample_id \
                 --localcores=$core \
                 --localmem=$mem \
                 --create-bam=true \
                 --nosecondary
