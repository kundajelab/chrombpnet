#!/bin/bash

cell_line=IMR90
data_type="DNASE_SE"
date=$(date +'%m.%d.%Y')
setting=$data_type"_"$date"_fold_0_hep_bias_transfer"
cur_file_name="imr90_fold_0retrain.sh"
### SIGNAL INPUT
fold=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/splits/fold_0.json
bias_h5=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/HEPG2/HEPG2_06.08.2022_bias_128_4_1234_0.8_fold_0/bias_model/bias.h5

overlap_peak=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR477RTP/preprocessing/downloads/peaks.bed.gz
chrom_sizes=reference/chrom.sizes
ref_fasta=reference/hg38.genome.fa

main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
neg_dir=$main_dir/negatives_data/
data_dir=$main_dir/data
output_dir=$main_dir/$setting
inputlen=2114
gpu=MIG-166d7783-762d-5f61-b31c-549eb4e0fba0


function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}


## CREATE DIRS
if [[ -d $main_dir ]] ; then
    echo "main director already exists"
else
    mkdir $main_dir
fi

if [[ -d $output_dir ]] ; then
    echo "output director already exists"
else
    mkdir $output_dir
fi



### STEP 2 - TRAIN CHROMBPNET MODEL

if [[ -d $output_dir/chrombpnet_model ]] ; then
    echo "skipping chrombpnet model training  - directory present "
else
    mkdir $output_dir/chrombpnet_model
    CUDA_VISIBLE_DEVICES=$gpu bash step6_train_chrombpnet_model.sh $ref_fasta $data_dir"/"$cell_line"_unstranded.bw" $overlap_peak $neg_dir/negatives_with_summit.bed $fold $bias_h5 $output_dir/chrombpnet_model $data_type
fi
