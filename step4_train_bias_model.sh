
reference_fasta=$1
bigwig_path=$2
overlap_peak=$3
nonpeaks=$4
fold=$5
bias_threshold_factor=$6
output_dir=$7

# defaults
inputlen=2114
outputlen=1000
filters=128
n_dilation_layers=4

# this script does the following -  
# (1) filters your peaks/nonpeaks (removes outliers and removes edge cases and creates a new filtered set)
# (2) filters non peaks based on the given bias threshold factor
# (3) Calculates the counts loss weight 
# (4) Creates a TSV file that can be loaded into the next step
python $PWD/src/helpers/hyperparameters/find_bias_hyperparams.py \
       --genome=$reference_fasta \
       --bigwig=$bigwig_path \
       --peaks=$overlap_peak \
       --nonpeaks=$nonpeaks \
       --outlier_threshold=0.99 \
       --chr_fold_path=$fold \
       --inputlen=$inputlen \
       --outputlen=$outputlen \
       --max_jitter=50 \
       --filters=$filters \
       --n_dilation_layers=$n_dilation_layers \
       --bias_threshold_factor=$bias_threshold_factor \
       --output_dir $output_dir 

# this script does the following -  
# (1) trains a model on the given peaks/nonpeaks
# (2) The parameters file input to this script should be TSV seperated 
python $PWD/src/training/train.py \
       --genome=$reference_fasta \
       --bigwig=$bigwig_path \
       --nonpeaks=$output_dir/filtered.bias_nonpeaks.bed \
       --params=$output_dir/bias_model_params.tsv \
       --output_prefix=$output_dir/bias \
       --chr_fold_path=$fold \
       --batch_size=64 \
       --architecture_from_file=$PWD/src/training/models/bpnet_model.py \
       --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss 

# predictions and metrics on the chrombpnet model trained
python $PWD/src/training/predict.py \
        --genome=$reference_fasta \
        --bigwig=$bigwig_path \
        --nonpeaks=$output_dir/filtered.bias_nonpeaks.bed \
        --chr_fold_path=$fold \
        --inputlen=$inputlen \
        --outputlen=$outputlen \
        --output_prefix=$output_dir/bias \
        --batch_size=256 \
        --model_h5=$output_dir/bias.h5 \
