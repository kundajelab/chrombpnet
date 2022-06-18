#!/bin/bash

cell_line=H1ESC
data_type="DNASE_SE"

date=$(date +'%m.%d.%Y')
cur_file_name="h1esc_dnase_se_deepshap.sh"
### SIGNAL INPUT

ref_fasta=reference/hg38.genome.fa

main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
data_dir=$main_dir/data
inputlen=2114
output_dir=$PWD/results/chrombpnet/$data_type/$cell_line/nautilus_runs_apr12/H1ESC_04.09.2022_bias_128_4_1234_0.8_fold_0/

gpu=0

function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}



if [[ -d $output_dir/chrombpnet_model/interpret ]] ; then
    echo "skipping bias model interpretation - directory present "
else
    mkdir $output_dir/chrombpnet_model/interpret
    logfile=$output_dir/chrombpnet_model/interpret/interpret.log
    touch $logfile

    echo $( timestamp ):"python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=results/chrombpnet/DNASE_SE/$cell_line/data/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/chrombpnet_model/interpret/$cell_line \
        --model_h5=$output_dir/chrombpnet_model/chrombpnet_wo_bias.h5" |  tee -a $logfile

    python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=results/chrombpnet/DNASE_SE/$cell_line/data/30K.subsample.overlap.bed \
        --output_prefix=$output_dir/chrombpnet_model/interpret/$cell_line \
        --model_h5=$output_dir/chrombpnet_model/chrombpnet_wo_bias.h5  | tee -a $logfile
fi
