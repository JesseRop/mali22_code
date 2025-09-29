#!/bin/bash
#BSUB -o /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir/%J.o
#BSUB -e /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir/%J.e
#BSUB -M 8000 -R "select[mem>8000] rusage[mem=8000]" -M 8000
#BSUB -q normal
#BSUB -n 2
 
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms8G -Xmx8G -Dnxf.pool.maxThreads=2000"
export NXF_VER=22.04.0-5697
  
# Load Nextflow module
module load nextflow/25.04.1-5946
  
# Run the Nextflow script with the specified parameters
nextflow run \
/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code_lnk/nf_scripts/submsn_intgrtn_pc.nf \
-w /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir/work \
-profile sanger \
-c /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir/profiles.config \
-qs 1000 \
-resume \
--softpack_module "HGI/softpack/groups/team222/Pf_scRNAseq/33" \
--scrpt "/nfs/users/nfs_j/jr35/multipurpose_scripts/seurat_sce/seu_integrtn_pc_auto.R" \
--input_seu "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/m2_all/v3_m2_asex_cr_filt_{all,noM55}_norm.RDS" \
--mthd_nm "harmony_PC_man" --cluster_resoln 0.4 --intgrtn_pc 8 --integrtd_pc_man "yes" \
--ncores "6" --mem "60 GB" --run_time="1.hour"  \
--findpc_mthd "perpendicular line" \
--mthd "HarmonyIntegration" \
--norm_method "LogNormalize" 


### 1. Run on merged seurat object after QC
# --input_seu "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/m2_all/all_msc_cb_no_dbt_mrg_noM12_sct.RDS" \
# --mthd_nm "harmony_PC_manl" --cluster_resoln 1 --intgrtn_pc 8 --integrtd_pc_man "yes" \
# --ncores "15" --mem "100 GB" --run_time="5.hour"  \

### 2. Run on QCd asexual
# --input_seu "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/m2_all/v3_m2_asex_cr_filt_{all,noM55}_norm.RDS" \
# --mthd_nm "harmony" --cluster_resoln 1 --intgrtn_pc 8 --integrtd_pc_man "no" \
# --ncores "6" --mem "60 GB" --run_time="1.hour"  \

### 2. Run on QCd asexual
# --input_seu "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/m2_all/v3_m2_asex_cr_filt_{all,noM55}_norm.RDS" \
# --mthd_nm "harmony_PC_man" --cluster_resoln 1 --intgrtn_pc 8 --integrtd_pc_man "yes" \
# --ncores "6" --mem "60 GB" --run_time="1.hour"  \

### 1. Run on asexuals after QC with PC specified manually
# --input_seu "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/m2_all/m2_lasse_asex_norm.RDS" \ 
# --mthd_nm "harmony" --cluster_resoln 0.4 --intgrtn_pc 8 --integrtd_pc_man "no" \
# --ncores "10" --mem "100 GB" --run_time="1.hour"

### 1. Run on asexuals after QC with PC specified manually
# --input_seu "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/m2_all/m2_lasse_asex_norm.RDS" \ 
# --mthd_nm "harmony_pc" --cluster_resoln 0.2 --intgrtn_pc 10 --integrtd_pc_man "yes" \
# --ncores "10" --mem "100 GB" --run_time="1.hour"

# ## clean up on exit 0 - delete this if you want to keep the work dir
# status=$?
# if [[ $status -eq 0 ]]; then
#   rm -r /path/to/some/dir/work
# fi
