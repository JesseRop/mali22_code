#!/bin/bash

##Submitting souporcell jobs

w_dir=/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data
# raw_dir=/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf_all_genes/
# raw_dir=/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf/
# raw_dir=/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/filtered_cellranger_h5/Pf/
# raw_dir_filt=/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/filtered_cellranger_h5/Pf_all_genes/
# raw_dir=/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs_rmHsapiens/Pf_all_genes/
raw_dir=/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf_all_genes/

## With or without human gene removal
# raw_dir=

## Job submission variables
# source_irods=("5736STDY11771536" "5736STDY11771545" "5736STDY11771544" "5736STDY11771535")
source_irods=( $(cut -d ',' -f2 $w_dir/raw/pf_solo_mixed_mdata_decode_cln.csv ) )
# source_irods=( $(cut -d ',' -f2 /nfs/users/nfs_j/jr35/Mali2/mali22_code/bash_scripts/pf_md2.csv))
# source_irods=("${source_irods[@]:1:2}")
# source_irods=("${source_irods[@]:28:1}")
# sample_name=( "msc3" "msc1272" "msc13" "msc1")
sample_name=( $(cut -d ',' -f1 $w_dir/raw/pf_solo_mixed_mdata_decode_cln.csv ) )
# sample_name=( $(cut -d ',' -f1 /nfs/users/nfs_j/jr35/Mali2/mali22_code/bash_scripts/pf_md2.csv))
# sample_name=("${sample_name[@]:1:2}")
# sample_name=("${sample_name[@]:28:1}")

##Expensive to rerun hence just runing for MSC14 and commenting out above code 
# source_irods=("5736STDY13782216")
# sample_name=( "MSC68")

# cb_sufx="dflt"
# cb_sufx="man_LRpt000005e250"
# cb_sufx="man_LRpt0001"
# cb_sufx="_rmHS_LRpt000005e200"

# l_rate=0.000005
l_rates=("0.0001" "0.000005")

# l_rate=0.00001
# epochs=200
epochs=("150" "250")

ncores=4
mem=80000

# exp_cells=6000
# exp_edrops=60000

for l_rate in "${l_rates[@]}"; do
    for epoch in "${epochs[@]}"; do
        for i in "${!sample_name[@]}";
        do
            echo "${sample_name[i]}"   
            scrpt=/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts_lnk/cbender_cellNo_edropsNo.sh
            input_raw_mtx=$raw_dir/${source_irods[i]}/outs/raw_feature_bc_matrix.h5
            # input_raw_mtx=$raw_dir_filt/${source_irods[i]}/outs/raw_feature_bc_matrix_rm_HS_gns.h5

            ##Original barcodes
            # o_dir=$w_dir/processed/Pf/${sample_name[i]}/cbender_${cb_sufx}/
            # o_dir=$w_dir/processed/Pf/${sample_name[i]}/cbender_custom/cb_LR${l_rate}_E${epoch}/
            o_dir=$w_dir/processed/Pf/${sample_name[i]}/cbender_custom_hs/cb_LR${l_rate}_E${epoch}/

            # o_file="${o_dir}/cb_${cb_sufx}.h5"
            o_file="${o_dir}/cb.h5"
            # exp_cells=$(cat "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/yascp_until_all/results/preprocessing/cellbender/${sample_name[i]}/cellbender-estimate_ncells_nemptydroplets/umi_count_estimates-expected_cells.txt")
            # exp_edrops=$(cat "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/yascp_until_all/results/preprocessing/cellbender/${sample_name[i]}/cellbender-estimate_ncells_nemptydroplets/umi_count_estimates-total_droplets_included.txt")

            exp_cells=$(cat "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/yascp_until_all_hs/results/preprocessing/cellbender/${sample_name[i]}/cellbender-estimate_ncells_nemptydroplets/umi_count_estimates-expected_cells.txt")
            exp_edrops=$(cat "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/yascp_until_all_hs/results/preprocessing/cellbender/${sample_name[i]}/cellbender-estimate_ncells_nemptydroplets/umi_count_estimates-total_droplets_included.txt")


            mkdir -p $o_dir/logs/
            
            eval $(echo "cd $o_dir & bash $scrpt $input_raw_mtx $o_dir $o_file $l_rate $epoch $exp_cells $exp_edrops $ncores $mem %J-%I")


        done
    done
done

