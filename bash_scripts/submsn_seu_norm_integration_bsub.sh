#!/bin/bash
#BSUB -o /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir/%J.o
#BSUB -e /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir/%J.e
#BSUB -M 8000 -R "select[mem>8000] rusage[mem=8000]" -M 8000
#BSUB -q normal
#BSUB -n 1
 
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms8G -Xmx8G -Dnxf.pool.maxThreads=2000"
export NXF_VER=22.04.0-5697
  
# Load Nextflow module
module load nextflow/25.04.1-5946
  
# Run the Nextflow script with the specified parameters
nextflow run \
/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code_lnk/nf_scripts/submsn_seu_norm.nf \
-w /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir/work \
-profile sanger \
-c /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir/profiles.config \
-qs 1000 \
-resume \
--softpack_module "HGI/softpack/groups/team222/Pf_scRNAseq/33" \
--scrpt "/nfs/users/nfs_j/jr35/multipurpose_scripts/seurat_sce/seu_stdnorm.R" \
--input_seu "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/m2_all/m2_lasse_ase.*.RDS" \
--hvg 400 --cluster_resoln 0.1 --cl_nm "seurat_clusters" --integrtd_pc30 "no" \
--ncores "8" --mem "80 GB" --run_time="5.min" \
--normlzn "norm" --umap_name "UMAP_" \
--umap_redctn_nm "umap.unintegrated" \
--redctn "pca" \
--pca_redctn_nm "pca" \
--resume


### 1. Run on merged seurat object after QC
# --input_seu "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/m2_all/all_msc_cb_no_dbt_mrg.RDS" \
# --hvg 700 --cluster_resoln 0.1 --cl_nm "seurat_clusters" --integrtd_pc30 "no" \
# --ncores "10" --mem "100 GB" --run_time="5.hour" \

### 1. Run on asexuals after QC with PC calculated by findPC
# --input_seu "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/m2_all/v3_m2_asex_cr_filt.RDS" \
# --hvg 500 --cluster_resoln 0.1 --cl_nm "seurat_clusters" --integrtd_pc30 "no" \
# --ncores "5" --mem "100 GB" --run_time="5.min" \

### 1. Run on asexuals after QC with PC calculated by findPC
# --input_seu "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/m2_all/m2_lasse_asex.RDS" \
# --hvg 500 --cluster_resoln 0.1 --cl_nm "seurat_clusters" --integrtd_pc30 "no" \
# --ncores "8" --mem "100 GB" --run_time="5.min" \

### 1. Run on asexuals after QC with PC specified manually
# --input_seu "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/m2_all/v3_m2_asex_cr_filt_{all,noM55}.RDS" \ 
# --hvg 400 --cluster_resoln 0.1 --cl_nm "seurat_clusters" --integrtd_pc30 "no" \
# --ncores "4" --mem "80 GB" --run_time="55.min" \

# --mthd_nm "harmony_pc_man" --cluster_resoln 0.4 --intgrtn_pc 6 --integrtd_pc_man "yes" \
# --ncores "10" --mem "100 GB" --run_time="1.hour"

### 1. Run on asexuals after QC with PC specified manually
# --input_seu "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf/m2_all/m2_lasse_ase.*.RDS" \
# --hvg 400 --cluster_resoln 0.1 --cl_nm "seurat_clusters" --integrtd_pc30 "no" \
# --ncores "8" --mem "80 GB" --run_time="5.min" \


# ## clean up on exit 0 - delete this if you want to keep the work dir
# status=$?
# if [[ $status -eq 0 ]]; then
#   rm -r /path/to/some/dir/work
# fi
