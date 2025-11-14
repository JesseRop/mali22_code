#!/bin/bash
#BSUB -o /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir/%J.o
#BSUB -e /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir/%J.e
#BSUB -M 8000 -R "select[mem>8000] rusage[mem=8000]" -M 8000
#BSUB -q oversubscribed
#BSUB -n 2
 
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms8G -Xmx8G -Dnxf.pool.maxThreads=2000"
export NXF_VER=22.04.0-5697
  
# Load Nextflow module
module load nextflow/25.04.1-5946
  
# Run the Nextflow script with the specified parameters
nextflow run \
/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code_lnk/souporcell/soupc_v25_m2.nf \
-w /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir/work_soupc \
-profile sanger \
-c /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir/profiles.config \
-qs 1000 \
-resume \
--id_decode "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/pf_all_id_decode.csv" \
--bam "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf_all_genes/5736STDY*/outs/possorted_genome_bam.bam" \
--bcodes "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/soupc_resc_hspf/barcodes.tsv.gz" \
--soup_dir "soupc_resc_hspf" \
--o_dir "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/" \
--scrpt "/nfs/users/nfs_j/jr35/multipurpose_scripts/soupc_v25_mnmap_hsat_310.sh" \
--ref_file "PlasmoDB-68_Pfalciparum3D7_Genome_VAR_masked.fasta" \
--hsat_ref_file "Pf68" \
--ncores "10" \
--mem "100 GB" \
--run_time "12.hour"
  
# ## clean up on exit 0 - delete this if you want to keep the work dir
# status=$?
# if [[ $status -eq 0 ]]; then
#   rm -r /path/to/some/dir/work
# fi

# 1. Run for Pf with only protein coding genes cellbender rescued cells
# --bam "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf/5736STDY*/outs/possorted_genome_bam.bam" \
# --bcodes "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/soup_resc/barcodes.tsv.gz" \
# --soup_dir "soupc_resc" \

# 1. Run for Pf with all genes in gtf cellbender rescued cells
# --bam "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf_all_genes/5736STDY*/outs/possorted_genome_bam.bam" \
# --bcodes "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/soupc_resc_hspf/barcodes.tsv.gz" \
# --soup_dir "soupc_resc_hspf" \
