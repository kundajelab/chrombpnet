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
logfile=$8
data_type=$7

# defaults
inputlen=2114
outputlen=1000
filters=512
n_dilation_layers=8
negative_sampling_ratio=0.1

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

echo $( timestamp ): "python $PWD/src/helpers/hyperparameters/find_hintatac_hyperparams.py \\
       --genome=$reference_fasta \\
       --bigwig=$bigwig_path \\
       --peaks=$overlap_peak \\
       --nonpeaks=$nonpeaks \\
       --outlier_threshold=0.99 \\
       --chr_fold_path=$fold \\
       --negative_sampling_ratio=$negative_sampling_ratio \\
       --inputlen=$inputlen \\
       --outputlen=$outputlen \\
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
echo $( timestamp ): "python $PWD/src/training/train.py \\
       --genome=$reference_fasta \\
       --bigwig=$bigwig_path \\
       --peaks=$output_dir/filtered.peaks.bed  \\
       --params=$output_dir/chrombpnet_model_params.tsv \\
       --output_prefix=$output_dir/hint_atac \\
       --chr_fold_path=$fold \\
       --batch_size=64 \\
       --architecture_from_file=$PWD/src/training/models/bpnet_model.py \\
       --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss" | tee -a $logfile

python $PWD/src/training/train.py \
       --genome=$reference_fasta \
       --bigwig=$bigwig_path \
       --peaks=$output_dir/filtered.peaks.bed  \
       --params=$output_dir/chrombpnet_model_params.tsv \
       --output_prefix=$output_dir/hint_atac \
       --chr_fold_path=$fold \
       --batch_size=64 \
       --architecture_from_file=$PWD/src/training/models/bpnet_model.py \
       --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss  | tee -a $logfile

# predictions and metrics on the bias model trained
echo $( timestamp ): "python $PWD/src/training/predict.py \\
        --genome=$reference_fasta \\
        --bigwig=$bigwig_path \\
        --peaks=$overlap_peak \\
        --chr_fold_path=$fold \\
        --inputlen=$inputlen \\
        --outputlen=$outputlen \\
        --output_prefix=$output_dir/hint_atac \\
        --batch_size=256 \\
        --model_h5=$output_dir/hint_atac.h5" | tee -a $logfile

python $PWD/src/training/predict.py \
        --genome=$reference_fasta \
        --bigwig=$bigwig_path \
        --peaks=$overlap_peak \
        --chr_fold_path=$fold \
        --inputlen=$inputlen \
        --outputlen=$outputlen \
        --output_prefix=$output_dir/hint_atac \
        --batch_size=256 \
        --model_h5=$output_dir/hint_atac.h5 | tee -a $logfile

# marginal footprinting
if [[ "$data_type" = "DNASE_SE" || "$data_type" = "DNASE_PE" ]] ; then
        echo $( timestamp ): "mkdir $output_dir/footprints" | tee -a $logfile
        mkdir $output_dir/footprints
        echo $( timestamp ): "python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \\
        -g $reference_fasta \\
        -r $output_dir/filtered.nonpeaks.bed \\
        --chr_fold_path=$fold \\
        -m $output_dir/hint_atac.h5 \\
        -bs 512 \\
        -o $output_dir/footprints/corrected \\
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \\
        -mo dnase_1,dnase_2" | tee -a $logfile
        python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \
        -g $reference_fasta \
        -r $output_dir/filtered.nonpeaks.bed \
        --chr_fold_path=$fold \
        -m $output_dir/hint_atac.h5 \
        -bs 512 \
        -o $output_dir/footprints/corrected \
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \
        -mo dnase_1,dnase_2 | tee -a $logfile
elif [[ "$data_type" = "ATAC_SE" || "$data_type" = "ATAC_PE"  ]] ; then
        echo $( timestamp ): "mkdir $output_dir/footprints" | tee -a $logfile
        mkdir $output_dir/footprints
        echo $( timestamp ): "python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \\
        -g $reference_fasta \\
        -r $output_dir/filtered.nonpeaks.bed \\
        --chr_fold_path=$fold \\
        -m $output_dir/hint_atac.h5 \\
        -bs 512 \\
        -o $output_dir/footprints/corrected \\
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \\
        -mo tn5_1,tn5_2,tn5_3,tn5_4,tn5_5" | tee -a $logfile
        python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \
        -g $reference_fasta \
        -r $output_dir/filtered.nonpeaks.bed \
        --chr_fold_path=$fold \
        -m $output_dir/hint_atac.h5 \
        -bs 512 \
        -o $output_dir/footprints/corrected \
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \
        -mo tn5_1,tn5_2,tn5_3,tn5_4,tn5_5 | tee -a $logfile
else
    echo "ERROR: unknown data type " $data_type | tee -a $logfile
fi

