#!/bin/bash

#BSUB -M 8000 -R "select[mem>8000] rusage[mem=8000]"
#BSUB -q oversubscribed
#BSUB -n 2
#BSUB -o /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir3/%J.edrops.o
#BSUB -e /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir3/%J.edrops.e
 
# Nextflow environment settings
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms8G -Xmx8G -Dnxf.pool.maxThreads=2000"
export NXF_VER=25.04.6-5954

# Ensure no conflicting nextflow module is loaded
module unload nextflow || true
module load nextflow/25.04.6-5954

# Base directories
BASE_DIR="/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2"
WORK_DIR="${BASE_DIR}/work_dir3"
DATA_DIR="${BASE_DIR}/data"
RAW_DIR="/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data"

# Optional parameter overrides (uncomment to customize)
LIMT="50,100,150,200"
RETAIN_NM="wo_retain,w_retain"
RETAIN_VAL="Inf,NULL"
FDR_L="0.01"

# Run the Nextflow edrops pipeline
nextflow run \
  ${BASE_DIR}/mali22_code_lnk/cell_calling/edrops.nf \
  -w ${WORK_DIR}/work_edrops \
  -profile sanger \
  -c ${WORK_DIR}/profiles.config \
  -qs 1000 \
  -resume \
  -with-dag edrops_dag.png \
  --id_decode "${DATA_DIR}/raw/pf_all_id_decode.csv" \
  --input_raw_mtx "${RAW_DIR}/{cellranger_runs,cellranger_runs_rmHsapiens}/Pf_all_genes/5736STDY*/outs/raw_feature_bc_matrix.h5" \
  --o_path "${DATA_DIR}/processed/Pf/" \
  --scrpt "/nfs/users/nfs_j/jr35/multipurpose_scripts/seurat_sce/edrops.R" \
  --bcrank_scrpt "/nfs/users/nfs_j/jr35/multipurpose_scripts/seurat_sce/bcrank_calc.R" \
  --qc_stats_scrpt "/nfs/users/nfs_j/jr35/multipurpose_scripts/seurat_sce/qc_stats.R" \
  --limt "${LIMT}" \
  --retain_nm "${RETAIN_NM}" \
  --retain_val "${RETAIN_VAL}" \
  --fdr_l "${FDR_L}"
