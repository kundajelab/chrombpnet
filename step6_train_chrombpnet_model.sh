reference_fasta=$1
bigwig_path=$2
overlap_peak=$3
nonpeaks=$4
fold=$5
bias_model=$6
output_dir=$7

# defaults
inputlen=2114
outputlen=1000
filters=512
n_dilation_layers=8
negative_sampling_ratio=0.1

# this script does the following -  
# (1) filters your peaks/nonpeaks (removes outliers and removes edge cases and creates a new filtered set)
# (2) scales the given bias model on the non-peaks
# (3) Calculates the counts loss weight 
# (4) Creates a TSV file that can be loaded into the next step
python $PWD/src/helpers/hyperparameters/find_chrombpnet_hyperparams.py \
       --genome=$reference_fasta \
       --bigwig=$bigwig_path \
       --peaks=$overlap_peak \
       --nonpeaks=$nonpeaks \
       --outlier_threshold=0.99 \
       --chr_fold_path=$fold \
       --inputlen=$inputlen \
       --outputlen=$outputlen \
       --max_jitter=10 \
       --filters=$filters \
       --n_dilation_layers=$n_dilation_layers \
       --bias_model_path=$bias_model \
       --output_dir=$output_dir 

# # this script does the following -  
# # (1) trains a model on the given peaks/nonpeaks
# # (2) The parameters file input to this script should be TSV seperatedp
python $PWD/src/training/train.py \
       --genome=$reference_fasta \
       --bigwig=$bigwig_path \
       --peaks=$output_dir/filtered.peaks.bed \
       --nonpeaks=$output_dir/filtered.nonpeaks.bed \
       --params=$output_dir/chrombpnet_model_params.tsv \
       --output_prefix=$output_dir/chrombpnet \
       --chr_fold_path=$fold \
       --epochs=40 \
       --batch_size=64 \
       --architecture_from_file=$PWD/src/training/models/chrombpnet_with_bias_model.py \
       --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss 

# # # predictions and metrics on the chrombpnet model trained
python $PWD/src/training/predict.py \
        --genome=$reference_fasta \
        --bigwig=$bigwig_path \
        --peaks=$output_dir/filtered.peaks.bed \
        --nonpeaks=$output_dir/filtered.nonpeaks.bed \
        --chr_fold_path=$fold \
        --inputlen=$inputlen \
        --outputlen=$outputlen \
        --output_prefix=$output_dir/chrombpnet \
        --batch_size=256 \
        --model_h5=$output_dir/chrombpnet.h5 \

# # # predictions and metrics on the chrombpnet model without bias trained
python $PWD/src/training/predict.py \
        --genome=$reference_fasta \
        --bigwig=$bigwig_path \
        --peaks=$output_dir/filtered.peaks.bed \
        --nonpeaks=$output_dir/filtered.nonpeaks.bed \
        --chr_fold_path=$fold \
        --inputlen=$inputlen \
        --outputlen=$outputlen \
        --output_prefix=$output_dir/chrombpnet_wo_bias \
        --batch_size=256 \
        --model_h5=$output_dir/chrombpnet_wo_bias.h5 \

# # # predictions and metrics on the bias model trained
python $PWD/src/training/predict.py \
        --genome=$reference_fasta \
        --bigwig=$bigwig_path \
        --peaks=$output_dir/filtered.peaks.bed \
        --nonpeaks=$output_dir/filtered.nonpeaks.bed \
        --chr_fold_path=$fold \
        --inputlen=$inputlen \
        --outputlen=$outputlen \
        --output_prefix=$output_dir/bias \
        --batch_size=256 \
        --model_h5=$output_dir/bias_model_scaled.h5 \

mkdir $output_dir/tn5_footprints
python src/evaluation/marginal_footprints/marginal_footprinting.py \
        -g $reference_fasta \
        -r $output_dir/filtered.nonpeaks.bed \
        -chr "chr1" \
        -m $output_dir/chrombpnet_wo_bias.h5 \
        -bs 512 \
        -o $output_dir/tn5_footprints/ \
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \
        -mo tn5_1,tn5_2,tn5_3,tn5_4,tn5_5

