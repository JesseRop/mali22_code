#!/bin/bash

##Submitting souporcell jobs

w_dir=/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data
raw_dir=/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf/

## Job submission variables
# source_irods=("5736STDY11771536" "5736STDY11771545" "5736STDY11771544" "5736STDY11771535")
source_irods=( $(cut -d ',' -f2 $w_dir/raw/pf_solo_mixed_mdata_decode.csv ) )
# source_irods=("${source_irods[@]:1:2}")
source_irods=("${source_irods[@]:28:1}")
# sample_name=( "msc3" "msc1272" "msc13" "msc1")
sample_name=( $(cut -d ',' -f1 $w_dir/raw/pf_solo_mixed_mdata_decode.csv ) )
# sample_name=("${sample_name[@]:1:2}")
sample_name=("${sample_name[@]:28:1}")

##Expensive to rerun hence just runing for MSC14 and commenting out above code 
# source_irods=("5736STDY11771546")
# sample_name=( "msc14")

ncores=2
mem=40000

for i in "${!sample_name[@]}";
do
    echo "${sample_name[i]}"   
    scrpt=/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts_lnk/cbender.sh
    input_raw_mtx=$raw_dir/${source_irods[i]}/outs/raw_feature_bc_matrix.h5

    ##Original barcodes
    o_dir=$w_dir/processed/Pf/${sample_name[i]}/cbender/
    o_file="${o_dir}/cb.h5"

    mkdir -p $o_dir/logs/
    
    eval $(echo "bash $scrpt $input_raw_mtx $o_dir $o_file $ncores $mem %J-%I")


done
