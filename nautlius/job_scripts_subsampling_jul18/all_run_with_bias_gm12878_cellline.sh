#!/bin/bash


date=$(date +'%m.%d.%Y')
bias_model=$1
seed=$2
cell_line=$3
fold_num=$4
data_type=$5

cur_file_name="all_run_with_bias.sh"

setting=$cell_line"_"$date"_bias_transfer_"$seed"_fold_"$fold_num"_data_type_"$data_type
fold=splits/"fold_"$fold_num".json"
negatives_bed=negatives_with_summit_$fold_num".bed"

### SIGNAL INPUT

if [[ "$fold_num" = "0"  ]] ; then
    negatives_bed=negatives_with_summit.bed
else
    negatives_bed=negatives_with_summit_$fold_num".bed"
fi

echo $negatives_bed

chrom_sizes=reference/chrom.sizes
ref_fasta=reference/hg38.genome.fa
output_dir=models/$setting
data_dir=$cell_line"_"$data_type/data
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




### STEP 2 - TRAIN CHROMBPNET MODEL

if [[ -d $output_dir/chrombpnet_model ]] ; then
    echo "skipping chrombpnet model training  - directory present "
else
    mkdir $output_dir/chrombpnet_model
    bash step6_train_chrombpnet_model.sh $ref_fasta $data_dir"/GM12878_unstranded.bw" $overlap_peak $data_dir/$negatives_bed $fold $bias_model $output_dir/chrombpnet_model $data_type $seed
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


