#!/bin/bash

cell_line=GM12878
data_type="ATAC_PE"

date=$(date +'%m.%d.%Y')
cur_file_name="gm12878_dil_layers_9_fold0.sh"
### SIGNAL INPUT

dil=9
setting=$data_type"_"$date"_dilation_layers_"$dil

blacklist_region=$PWD/reference/GRch38_unified_blacklist.bed.gz
chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
fold=$PWD/splits/fold_0.json

oak_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/$data_type/$cell_line/

main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
data_dir=$main_dir/data
output_dir=$main_dir/$setting
neg_dir=$main_dir/negatives_data

inputlen=4114
gpu=0

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



overlap_peak=$data_dir/peaks_no_blacklist.bed


transfer_bias_model=$PWD/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/bias_model/bias_input_4114.h5

### STEP 2 - TRAIN CHROMBPNET MODEL

if [[ -d $output_dir/chrombpnet_model ]] ; then
    echo "skipping chrombpnet model training  - directory present "
else
    mkdir $output_dir/chrombpnet_model
    CUDA_VISIBLE_DEVICES=$gpu bash step6_train_chrombpnet_model.sh $ref_fasta $data_dir"/"$cell_line"_unstranded.bw" $overlap_peak $neg_dir/negatives_with_summit.bed $fold $transfer_bias_model $output_dir/chrombpnet_model $data_type 1234 $dil
fi

