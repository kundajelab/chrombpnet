# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

reference_fasta=$1
bigwig_path=$2
overlap_peak=$3
nonpeaks=$4
fold=$5
output_dir=$6
data_type=$7
bias_bigwig=$8
logfile=$9


bias_filters=128
bias_dil=4
seed=1234

# defaults
inputlen=2114
outputlen=1000
bias_threshold_factor=0.4
filters=$bias_filters
n_dilation_layers=$bias_dil
seed=$seed

if [ -z "$bias_filters" ]
  then
    filters=128
fi

if [ -z "$bias_dil" ]
  then
    n_dilation_layers=4
fi

if [ -z "$seed" ]
  then
    seed=1234
fi

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

negative_sampling_ratio=0.1

echo $( timestamp ): "python $PWD/src/helpers/hyperparameters/find_hintatac_hyperparams.py \\
       --genome=$reference_fasta \\
       --bigwig=$bigwig_path \\
       --peaks=$overlap_peak \\
       --nonpeaks=$nonpeaks \\
       --outlier_threshold=0.99 \\
       --chr_fold_path=$fold \\
       --inputlen=$inputlen \\
       --outputlen=$outputlen \\
       --negative_sampling_ratio=$negative_sampling_ratio \\
       --max_jitter=500 \\
       --filters=$filters \\
       --n_dilation_layers=$n_dilation_layers \\
       --output_dir=$output_dir " | tee -a $logfile

python $PWD/src/helpers/hyperparameters/find_hintatac_hyperparams.py \
       --genome=$reference_fasta \
       --bigwig=$bigwig_path \
       --peaks=$overlap_peak \
       --nonpeaks=$nonpeaks \
       --outlier_threshold=0.99 \
       --chr_fold_path=$fold \
       --negative_sampling_ratio=$negative_sampling_ratio \
       --inputlen=$inputlen \
       --outputlen=$outputlen \
       --max_jitter=500 \
       --filters=$filters \
       --n_dilation_layers=$n_dilation_layers \
       --output_dir=$output_dir | tee -a $logfile

# this script does the following -  
# (1) trains a model on the given peaks/nonpeaks
# (2) The parameters file input to this script should be TSV seperated 
echo $( timestamp ): "python $PWD/src/training/with_bigwig_train.py \\
       --genome=$reference_fasta \\
       --bigwig=$bigwig_path \\
       --bias_bigwig=$bias_bigwig \\
       --peaks=$output_dir/filtered.peaks.bed \\
       --params=$output_dir/chrombpnet_model_params.tsv \\
       --output_prefix=$output_dir/tobias_bias \\
       --chr_fold_path=$fold \\
       --seed=$seed \\
       --batch_size=64 \\
       --architecture_from_file=$PWD/src/training/models/smooth_bigwig.py \\
       --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss" | tee -a $logfile

python $PWD/src/training/with_bigwig_train.py \
       --genome=$reference_fasta \
       --bigwig=$bigwig_path \
       --bias_bigwig=$bias_bigwig \
       --peaks=$output_dir/filtered.peaks.bed \
       --params=$output_dir/chrombpnet_model_params.tsv \
       --output_prefix=$output_dir/tobias_bias \
       --chr_fold_path=$fold \
       --seed=$seed \
       --batch_size=64 \
       --architecture_from_file=$PWD/src/training/models/smooth_bigwig.py \
       --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss  | tee -a $logfile

