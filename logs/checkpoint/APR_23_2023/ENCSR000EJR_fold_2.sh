#!/bin/bash

cell_line=ENCSR000EJR
data_type="DNASE_SE"

date=$(date +'%m.%d.%Y')
#setting=$data_type"_"$date"_invivo_bias"
cur_file_name="ENCSR000EJR_fold_2.sh"
setting=DNASE_SE_03.16.2023_1.0_fold_2

### SIGNAL INPUT

in_bigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR000EJR/preprocessing/bigWigs/ENCSR000EJR.bigWig
overlap_peak=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR000EJR/preprocessing/downloads/peaks.bed.gz
chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
blacklist_region=/mnt/data/annotations/blacklist/GRch38/GRch38_unified_blacklist.bed.gz
neg_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/ENCSR000EJR/negatives_data_2/
#transfer_bias_model=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/ENCSR000EKP/DNASE_SE_02.13.2023_0.9/bias_model/bias.h5
#transfer_bias_model=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/ENCSR000EKP/DNASE_SE_02.17.2023_1.0/bias_model/bias.h5
#transfer_bias_model=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/ENCSR000EJR/DNASE_SE_03.13.2023_0.8/bias_model/bias.h5
transfer_bias_model=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/ENCSR000EJR/DNASE_SE_03.16.2023_1.0_fold_2/bias_model/bias.h5
fold=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/splits/fold_2.json
gpu=2

oak_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/$data_type/$cell_line/
mkdir $oak_dir
main_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/$cell_line/
data_dir=$main_dir/data/

output_dir=$main_dir/$setting

inputlen=2114


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
    CUDA_VISIBLE_DEVICES=$gpu bash step6_train_chrombpnet_model.sh $ref_fasta $in_bigwig $overlap_peak $neg_dir/negatives_with_summit.bed $fold $transfer_bias_model $output_dir/chrombpnet_model $data_type
fi


#zcat $overlap_peak | shuf -n 30000 > $data_dir/30K.subsample.overlap.bed

### INTERPRET SEQUENCE MODEL

if [[ -f $output_dir/chrombpnet_model/interpret/$cell_line.counts_scores.h5 ]] ; then
    echo "skipping sequence model interpretation - directory present "
else
    mkdir $output_dir/chrombpnet_model/interpret

    logfile=$output_dir/chrombpnet_model/interpret/interpret.log
    touch $logfile

    echo $( timestamp ):"python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=$data_dir/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/chrombpnet_model/interpret/$cell_line \
        --model_h5=$output_dir/chrombpnet_model/chrombpnet_wo_bias.h5" |  tee -a $logfile

    CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=$data_dir/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/chrombpnet_model/interpret/$cell_line \
        --model_h5=$output_dir/chrombpnet_model/chrombpnet_wo_bias.h5  | tee -a $logfile
fi


if [[ -d $oak_dir/$setting/ ]]; then
    echo "dir exists"
else
    mkdir $oak_dir/$setting/
    mkdir $oak_dir/$setting/SIGNAL
fi



if [[ -f $oak_dir/$setting/SIGNAL/$cell_line.counts_scores.h5  && -f $oak_dir/$setting/SIGNAL/$cell_line.profile_scores.h5 ]] ; then
    echo "chrombpnet model files already present on oak"
else
    echo "copying counts and profile scores to oak"
    cp $output_dir/chrombpnet_model/interpret/$cell_line.counts_scores.h5 $oak_dir/$setting/SIGNAL/
    cp $output_dir/chrombpnet_model/interpret/$cell_line.profile_scores.h5 $oak_dir/$setting/SIGNAL/
fi

