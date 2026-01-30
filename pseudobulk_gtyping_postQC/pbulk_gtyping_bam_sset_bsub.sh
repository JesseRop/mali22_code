#!/bin/bash
#BSUB -o /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir3/%J.o
#BSUB -e /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir3/%J.e
#BSUB -M 8000 -R "select[mem>8000] rusage[mem=8000]" -M 8000
#BSUB -q oversubscribed
#BSUB -n 2
 
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms8G -Xmx8G -Dnxf.pool.maxThreads=2000"
export NXF_VER=22.04.0-5697
  
# Load Nextflow module
module load nextflow/25.04.6-5954
  
# Run the Nextflow script with the specified parameters
nextflow run \
/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code_lnk/nf_scripts/pbulk_gtyping_bam_sset.nf \
-w /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir3/work \
-profile sanger \
-c /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir3/profiles.config \
-qs 1000 \
-resume \
--sv_dir "pbulk_gtypes_cln_1VCF3" \
--bam "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/soupc_resc_hspf3/*/parent/out_dir/souporcell_minimap_tagged_sorted.bam" \
--bcodes "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/pbulk_gtypes_cln_1VCF3/bcodes/minmap/st*/*tsv" \
--odir "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/" \
--mapper "minmap" \
--scrpt "/nfs/users/nfs_j/jr35/multipurpose_scripts/cr_subset_bam_linux.sh" \
--save_bam_sset true 


### 1. Quick run on MSC68 and MSC70 with downsampling and dry run of vartrix
# --sv_dir "pbulk_gtypes_cln_1VCF" \
# --bam "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/{MSC68,MSC70}/pbulk_gtypes_cln_1VCF/str*/new_bams/*SC*_sset_sorted_rg_don.bam" \
# --bcodes "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/pbulk_gtypes_cln_1VCF/bcodes/minmap/str*/*tsv" \
# --downsample_bam "true" \
# --vartrix_dryrun "true"

### 1. On Jan28 2026 for everything post combining jumpcode
# --sv_dir "pbulk_gtypes_cln_1VCF3" \
# --bam "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/{MSC68,MSC70}/pbulk_gtypes_cln_1VCF3/str*/new_bams/*SC*_sset_sorted_rg_don.bam" \
# --bcodes "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/MSC*/pbulk_gtypes_cln_1VCF3/bcodes/minmap/str*/*tsv" \
# --downsample_bam "true" \
# --vartrix_dryrun "true"

# ## clean up on exit 0 - delete this if you want to keep the work dir
# status=$?
# if [[ $status -eq 0 ]]; then
#   rm -r /path/to/some/dir/work
# fi
