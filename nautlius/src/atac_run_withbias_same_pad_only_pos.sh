#!/bin/bash


date=$(date +'%m.%d.%Y')
seed=$1
cell_line=$2
fold_num=$3
dillayers=$4
bias_model_name=$cell_line/$5
inputlen=2114

echo $bias_model_name
data_type="ATAC_PE"
cur_file_name="atac_run_withbias_same_pad_only_pos.sh"
setting=$cell_line"_"$date"_"$seed"_"$dillayers"_"$inputlen"_"$fold_num"_samepad_onlypos"
fold=splits/"fold_"$fold_num".json"

### SIGNAL INPUT


chrom_sizes=reference/chrom.sizes
ref_fasta=reference/hg38.genome.fa
output_dir=models/$setting
data_dir=$cell_line/data
overlap_peak=$data_dir/peaks_no_blacklist.bed

if [[ "$fold_num" = "0"  ]] ; then
    negatives_bed=negatives_with_summit.bed
else
    negatives_bed=negatives_with_summit_$fold_num".bed"
fi

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
    bash step6_train_chrombpnet_model_dillayer_same_pad_only_pos.sh $ref_fasta $data_dir"/"$cell_line"_unstranded.bw" $overlap_peak $data_dir/$negatives_bed $fold $bias_model_name $output_dir/chrombpnet_model $data_type $seed $dillayers $inputlen
fi

cp -r /scratch/chrombpnet/models/* /chrombpnet/ATAC_PE/$cell_line/

### INTERPRET SEQUENCE MODEL

#if [[ -d $output_dir/chrombpnet_model/interpret ]] ; then
#    echo "skipping sequence model interpretation - directory present "
#else
#    mkdir $output_dir/chrombpnet_model/interpret

#    logfile=$output_dir/chrombpnet_model/interpret/interpret.log
#    touch $logfile


#    echo $( timestamp ):"python $PWD/src/evaluation/interpret/interpret.py \
#        --genome=$ref_fasta \
#        --regions=$data_dir/30K.subsample.overlap.bed \
#        --output_prefix=$output_dir/chrombpnet_model/interpret/$cell_line \
#        --model_h5=$output_dir/chrombpnet_model/chrombpnet_wo_bias.h5" |  tee -a $logfile

#    python $PWD/src/evaluation/interpret/interpret.py \
#        --genome=$ref_fasta \
#        --regions=$data_dir/30K.subsample.overlap.bed \
#        --output_prefix=$output_dir/chrombpnet_model/interpret/$cell_line \
#        --model_h5=$output_dir/chrombpnet_model/chrombpnet_wo_bias.h5  | tee -a $logfile
#fi

#cp -r /scratch/chrombpnet/models/* /chrombpnet/ATAC_PE/$cell_line/


