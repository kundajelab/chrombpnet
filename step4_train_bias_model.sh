#!/bin/bash
# exit when any command fails
set -e
set -o pipefail

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command failed with exit code $?."' EXIT

reference_fasta=${1?param missing - reference_fasta}
bigwig_path=${2?param missing - bigwig_path}
overlap_peak=${3?param missing - overlap_peak}
nonpeaks=${4?param missing - nonpeaks}
fold=${5?param missing - fold}
bias_threshold_factor=${6?param missing - bias_threshold_factor}
output_dir=${7?param missing - output_dir}
filters=${8:-128}
n_dilation_layers=${9?:-4}
seed=${10:-1234}
logfile=$11

# defaults
inputlen=2114
outputlen=1000

function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

# create the log file
if [ -z "$logfile" ]
  then
    echo "No logfile supplied - creating one"
    logfile=$output_dir"/train_bias_model.log"
    touch $logfile
fi

# this script does the following -  
# (1) filters your peaks/nonpeaks (removes outliers and removes edge cases and creates a new filtered set)
# (2) filters non peaks based on the given bias threshold factor
# (3) Calculates the counts loss weight 
# (4) Creates a TSV file that can be loaded into the next step
echo $( timestamp ): "chrombpnet_bias_hyperparams \\
       --genome=$reference_fasta \\
       --bigwig=$bigwig_path \\
       --peaks=$overlap_peak \\
       --nonpeaks=$nonpeaks \\
       --outlier_threshold=0.99 \\
       --chr_fold_path=$fold \\
       --inputlen=$inputlen \\
       --outputlen=$outputlen \\
       --max_jitter=0 \\
       --filters=$filters \\
       --n_dilation_layers=$n_dilation_layers \\
       --bias_threshold_factor=$bias_threshold_factor \\
       --output_dir $output_dir" | tee -a $logfile

chrombpnet_bias_hyperparams \
    --genome=$reference_fasta \
    --bigwig=$bigwig_path \
    --peaks=$overlap_peak \
    --nonpeaks=$nonpeaks \
    --outlier_threshold=0.99 \
    --chr_fold_path=$fold \
    --inputlen=$inputlen \
    --outputlen=$outputlen \
    --max_jitter=0 \
    --filters=$filters \
    --n_dilation_layers=$n_dilation_layers \
    --bias_threshold_factor=$bias_threshold_factor \
    --output_dir $output_dir | tee -a $logfile

# this script does the following -  
# (1) trains a model on the given peaks/nonpeaks
# (2) The parameters file input to this script should be TSV seperated
bpnet_model_path=`which bpnet_model.py`
echo $( timestamp ): "chrombpnet_train \\
       --genome=$reference_fasta \\
       --bigwig=$bigwig_path \\
       --nonpeaks=$output_dir/filtered.bias_nonpeaks.bed \\
       --params=$output_dir/bias_model_params.tsv \\
       --output_prefix=$output_dir/bias \\
       --chr_fold_path=$fold \\
       --seed=$seed \\
       --batch_size=64 \\
       --architecture_from_file=$bpnet_model_path \\
       --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss" | tee -a $logfile

chrombpnet_train \
    --genome=$reference_fasta \
    --bigwig=$bigwig_path \
    --nonpeaks=$output_dir/filtered.bias_nonpeaks.bed \
    --params=$output_dir/bias_model_params.tsv \
    --output_prefix=$output_dir/bias \
    --chr_fold_path=$fold \
    --seed=$seed \
    --batch_size=64 \
    --architecture_from_file=$bpnet_model_path \
    --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss  | tee -a $logfile

# predictions and metrics on the bias model trained
echo $( timestamp ): "chrombpnet_predict \\
        --genome=$reference_fasta \\	 
        --bigwig=$bigwig_path \\  
        --nonpeaks=$output_dir/filtered.bias_nonpeaks.bed \\
        --chr_fold_path=$fold \\
        --inputlen=$inputlen \\
        --outputlen=$outputlen \\
        --output_prefix=$output_dir/bias \\
        --batch_size=256 \\
        --model_h5=$output_dir/bias.h5" | tee -a $logfile

chrombpnet_predict \
    --genome=$reference_fasta \
    --bigwig=$bigwig_path \
    --nonpeaks=$output_dir/filtered.bias_nonpeaks.bed \
    --chr_fold_path=$fold \
    --inputlen=$inputlen \
    --outputlen=$outputlen \
    --output_prefix=$output_dir/bias \
    --batch_size=256 \
    --model_h5=$output_dir/bias.h5 | tee -a $logfile
