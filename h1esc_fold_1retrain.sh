#!/bin/bash

cell_line=H1ESC
data_type="DNASE_SE"
date=$(date +'%m.%d.%Y')
setting=$data_type"_"$date"_fold_1"
cur_file_name="h1esc_fold_1retrain.sh"
### SIGNAL INPUT
fold=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/splits/fold_1.json
bias_h5=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/H1ESC/H1ESC_07.07.2022_bias_128_4_1234_0.8_fold_1_data_type_DNASE_SE/bias_model/bias.h5

overlap_peak=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR000EMU/preprocessing/downloads/peaks.bed.gz
blacklist_region=/mnt/data/annotations/blacklist/GRch38/GRch38_unified_blacklist.bed.gz
chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
neg_dir=$main_dir/negatives_data_1/
data_dir=$main_dir/data
output_dir=$main_dir/$setting
inputlen=2114
gpu=1


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
