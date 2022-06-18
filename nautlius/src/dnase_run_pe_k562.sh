#!/bin/bash


date=$(date +'%m.%d.%Y')
bias_filters=$1
bias_dil=$2
seed=$3
bias_threshold_factor=$4
cell_line=$5
fold_num=$6

data_type="DNASE_PE"
cur_file_name="dnase_run_pe_k562.sh"
setting=$cell_line"_"$date"_bias_"$bias_filters"_"$bias_dil"_"$seed"_"$bias_threshold_factor"_"$fold_num
fold=splits/$fold_num".json"

### SIGNAL INPUT


chrom_sizes=reference/chrom.sizes
ref_fasta=reference/hg38.genome.fa
output_dir=models/$setting
data_dir=$cell_line/data
overlap_peak=$data_dir/peaks_no_blacklist.bed
inputlen=2114


function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}


if [[ -d $output_dir ]] ; then
    echo "output director already exists"
else
    mkdir models
    mkdir $output_dir
fi



### STEP 1 - TRAIN BIAS MODEL
if [[ -d $output_dir/bias_model ]] ; then
    echo "skipping bias model training  - directory present "
else
    mkdir $output_dir/bias_model
    bash step4_train_bias_model.sh $ref_fasta $data_dir"/"$cell_line"_unstranded.bw" $overlap_peak $data_dir/negatives_with_summit.bed $fold $bias_threshold_factor $output_dir/bias_model $bias_filters $bias_dil $seed
fi

cp -r /scratch/chrombpnet/models/* /chrombpnet/$data_type/$cell_line/

### STEP 2 - TRAIN CHROMBPNET MODEL

if [[ -d $output_dir/chrombpnet_model ]] ; then
    echo "skipping chrombpnet model training  - directory present "
else
    mkdir $output_dir/chrombpnet_model
    bash step6_train_chrombpnet_model.sh $ref_fasta $data_dir"/"$cell_line"_unstranded.bw" $overlap_peak $data_dir/negatives_with_summit.bed $fold $output_dir/bias_model/bias.h5 $output_dir/chrombpnet_model $data_type $seed
fi

cp -r /scratch/chrombpnet/models/* /chrombpnet/$data_type/$cell_line/

### INTERPRET SEQUENCE MODEL

if [[ -d $output_dir/chrombpnet_model/interpret ]] ; then
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

    python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=$data_dir/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/chrombpnet_model/interpret/$cell_line \
        --model_h5=$output_dir/chrombpnet_model/chrombpnet_wo_bias.h5  | tee -a $logfile
fi

cp -r /scratch/chrombpnet/models/* /chrombpnet/$data_type/$cell_line/

### INTERPRET BIAS MODEL


if [[ -d $output_dir/bias_model/interpret ]] ; then
    echo "skipping bias model interpretation - directory present "
else
    mkdir $output_dir/bias_model/interpret
    logfile=$output_dir/bias_model/interpret/interpret.log
    touch $logfile

    echo $( timestamp ):"python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=$data_dir/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/bias_model/interpret/$cell_line \
        --model_h5=$output_dir/bias_model/bias.h5" |  tee -a $logfile

    python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=$data_dir/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/bias_model/interpret/$cell_line \
        --model_h5=$output_dir/bias_model/bias.h5  | tee -a $logfile
fi

cp -r /scratch/chrombpnet/models/* /chrombpnet/$data_type/$cell_line/

