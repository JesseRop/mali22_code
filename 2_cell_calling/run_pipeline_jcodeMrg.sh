#!/bin/bash

#BSUB -M 8000 -R "select[mem>8000] rusage[mem=8000]"
#BSUB -q oversubscribed
#BSUB -n 2
#BSUB -o /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir2/%J.pipeline.o
#BSUB -e /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/work_dir2/%J.pipeline.e
 
# Nextflow environment settings
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms8G -Xmx8G -Dnxf.pool.maxThreads=2000"
export NXF_VER=25.04.6-5954

# Ensure no conflicting nextflow module is loaded
module unload nextflow || true
module load nextflow/25.04.6-5954

# Base directories
BASE_DIR="/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2"
WORK_DIR="${BASE_DIR}/work_dir2"
DATA_DIR="${BASE_DIR}/data"
RAW_DIR="/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data"

# Set parameters
L_RATES="0.0001,0.000005"
EPOCHS="150,250"
SPECIES_CRANGER="cellranger_runs_rmHsapiens cellranger_runs"
SPECIES_CBENDER="_rmHs _wHsPf"

# Run the Nextflow pipeline
nextflow run \
  ${BASE_DIR}/mali22_code_lnk/cell_calling/main.nf \
  -w ${WORK_DIR}/work_pipeline \
  -profile sanger \
  -c ${WORK_DIR}/profiles.config \
  -qs 1000 \
  -resume \
  -with-dag pipeline.png \
  --outdir "${DATA_DIR}" \
  --input_tenx "${RAW_DIR}/{cellranger_runs,cellranger_runs_rmHsapiens}/Pf_all_genes/5736STDY*/outs/raw_feature_bc_matrix" \
  --decode_csv "${DATA_DIR}/raw/pf_solo_mixed_jcode_merge_decode.csv" \
  --l_rates "${L_RATES}" \
  --epochs "${EPOCHS}" \
  \
  --script_umi "${BASE_DIR}/mali22_code_lnk/cell_calling/get_estimates_from_umi_counts.py" \
  --cellbender_script "${BASE_DIR}/../multipurpose_scripts_lnk/cbender_cellNo_edropsNo.sh" \
  --metadata_script "${BASE_DIR}/mali22_code_lnk/python_nbooks/cbender_mdata_csv.py" \
  --yascp_utils "/software/team222/jr35/yascp_utils" \
  --cellbender_env "/software/team222/jr35/cellbender_gpu/cellbender_gpu_env" \
  \
  --estimate_params_umis.expected_nemptydroplets_umi_cutoff "0" \
  --estimate_params_umis.method_estimate_ncells "dropletutils::barcoderanks::inflection" \
  --estimate_params_umis.lower_bound_umis_estimate_ncells "251" \
  --estimate_params_umis.method_estimate_nemptydroplets "dropletutils::barcoderanks::inflection,dropletutils::barcoderanks::knee,0.33" \
  --estimate_params_umis.lower_bound_umis_estimate_nemptydroplets "10" \
  --estimate_params_umis.upper_bound_umis_estimate_nemptydroplets "250" \
  --estimate_params_umis.estimate_nemptydroplets_umi_add_factor "0" \
  --estimate_params_umis.estimate_nemptydroplets_subtract_cell_factor "0" \
  --estimate_params_umis.min_droplets "0" \
  \
  --ncores "3" \
  --mem "30000" \
  --species_cranger ${SPECIES_CRANGER} \
  --species_cbender ${SPECIES_CBENDER} \
  \
  --obs_fname "cb_mdata.csv" \
  --var_fname "cb_gene_mdata.csv" #\
#  --test_sample "5736STDY11771536" ## To only run a single sample, MSC3 for testing or for running select samples for various reasons

# Exit status handling
# status=$?
# if [[ $status -ne 0 ]]; then
#   echo "Pipeline failed with exit code $status"
#   exit $status
# fi

## Dry run
#   \
#  --test_sample "5736STDY11771536"

# L_RATES="0.0001"
# EPOCHS="20"
# SPECIES_CRANGER="cellranger_runs_rmHsapiens"
# SPECIES_CBENDER="_rmHs"

## Default run

# L_RATES="0.0001,0.000005"
# EPOCHS="150,250"
# SPECIES_CRANGER="cellranger_runs_rmHsapiens cellranger_runs"
# SPECIES_CBENDER="_rmHs _wHsPf"