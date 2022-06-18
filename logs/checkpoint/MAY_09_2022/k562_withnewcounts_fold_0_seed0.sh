#!/bin/bash

cell_line=K562
data_type="ATAC_PE"

cur_file_name="k562_withnewcounts_fold_0_seed0.sh"
setting=K562_02.08.2022_bias_128_4_1234
### SIGNAL INPUT


blacklist_region=$PWD/reference/GRch38_unified_blacklist.bed.gz
chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
fold=$PWD/splits/fold_0.json

oak_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/$data_type/$cell_line/

main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
data_dir=$main_dir/data
output_dir=$PWD/results/chrombpnet/ATAC_PE/K562_stability/K562_02.08.2022_bias_128_4_1234/
neg_dir=$PWD/results/chrombpnet_feb_04/ATAC_PE/K562/negatives_data

inputlen=2114
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


transfer_bias_model=$PWD/results/chrombpnet/ATAC_PE/K562_stability/K562_02.08.2022_bias_128_4_1234/bias_model/bias.h5

### STEP 2 - TRAIN CHROMBPNET MODEL

if [[ -d $output_dir/chrombpnet_model_with_new_counts ]] ; then
    echo "skipping chrombpnet model training  - directory present "
else
    mkdir $output_dir/chrombpnet_model_with_new_counts
    CUDA_VISIBLE_DEVICES=$gpu bash step6_train_chrombpnet_model_newcounts.sh $ref_fasta $data_dir"/"$cell_line"_unstranded.bw" $overlap_peak $neg_dir/negatives_with_summit.bed $fold $transfer_bias_model $output_dir/chrombpnet_model_with_new_counts $data_type
fi

### INTERPRET SEQUENCE MODEL

if [[ -d $output_dir/chrombpnet_model_with_new_counts/interpret ]] ; then
    echo "skipping sequence model interpretation - directory present "
else
    mkdir $output_dir/chrombpnet_model_with_new_counts/interpret

    logfile=$output_dir/chrombpnet_model_with_new_counts/interpret/interpret.log
    touch $logfile

    echo $( timestamp ):"python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=$data_dir/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/chrombpnet_model_with_new_counts/interpret/$cell_line \
        --model_h5=$output_dir/chrombpnet_model_with_new_counts/chrombpnet_wo_bias.h5" |  tee -a $logfile

    CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=$data_dir/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/chrombpnet_model_with_new_counts/interpret/$cell_line \
        --model_h5=$output_dir/chrombpnet_model_with_new_counts/chrombpnet_wo_bias.h5  | tee -a $logfile
fi

if [[ -d $oak_dir/$setting/ ]]; then
    echo "dir exists"
else
    mkdir $oak_dir/$setting/
fi

mkdir $oak_dir/$setting/SIGNAL_newcounts

if [[ -f $oak_dir/$setting/SIGNAL_newcounts/$cell_line.counts_scores.h5  && -f $oak_dir/$setting/SIGNAL_newcounts/$cell_line.profile_scores.h5 ]] ; then
    echo "chrombpnet model files already present on oak"
else
    echo "copying counts and profile scores to oak"
    cp $output_dir/chrombpnet_model_with_new_counts/interpret/$cell_line.counts_scores.h5 $oak_dir/$setting/SIGNAL_newcounts/
    cp $output_dir/chrombpnet_model_with_new_counts/interpret/$cell_line.profile_scores.h5 $oak_dir/$setting/SIGNAL_newcounts/
fi

