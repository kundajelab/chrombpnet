#!/bin/bash


date=$(date +'%m.%d.%Y')
seed=1234
cell_line=K562
fold_num=0
dillayers=8
#inputlen=8500
inputlen=2114
gpu=0

bias_mod=HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0
bias_model_name=results/chrombpnet/ATAC_PE/HEPG2/nautilus_runs_may20/$bias_mod/bias_model/bias.h5

echo $bias_model_name
data_type="ATAC_PE"
cur_file_name="atac_run_withbias_only_pos_brahma.sh"
setting=$cell_line"_"$date"_"$seed"_"$dillayers"_"$inputlen"_"$fold_num"_only_pos"
fold=splits/"fold_"$fold_num".json"
### SIGNAL INPUT


chrom_sizes=reference/chrom.sizes
ref_fasta=reference/hg38.genome.fa
output_dir=results/chrombpnet/ATAC_PE/$cell_line/$setting
data_dir=results/chrombpnet/ATAC_PE/$cell_line/data
overlap_peak=$data_dir/peaks_no_blacklist.bed

if [[ "$fold_num" = "0"  ]] ; then
    negatives_bed=results/chrombpnet/ATAC_PE/HEPG2/negatives_data/negatives_with_summit.bed
else
    negatives_bed=results/chrombpnet/ATAC_PE/HEPG2/negatives_data_$fold_num/negatives_with_summit".bed"
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
    mkdir $output_dir
fi



### STEP 2 - TRAIN CHROMBPNET MODEL

if [[ -d $output_dir/chrombpnet_model ]] ; then
    echo "skipping chrombpnet model training  - directory present "
else
    mkdir $output_dir/chrombpnet_model
    CUDA_VISIBLE_DEVICES=$gpu bash step6_train_chrombpnet_model_dillayer_same_pad_only_pos.sh $ref_fasta $data_dir"/"$cell_line"_unstranded.bw" $overlap_peak $negatives_bed $fold $bias_model_name $output_dir/chrombpnet_model $data_type $seed $dillayers $inputlen
fi


#cp -r /scratch/chrombpnet/models/* /chrombpnet/ATAC_PE/$cell_line/

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


