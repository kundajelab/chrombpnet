
reference_fasta=$1
bigwig_path=$2
overlap_peak=$3
nonpeaks=$4
fold=$5
bias_model=$6

# defaults
inputlen=2114
outputlen=1000
filters=512
n_dilation_layers=8
negative_sampling_ratio=0.1

mkdir output/chrombpnet_model
output_dir=output/chrombpnet_model

#creates a bias_params.txt file in the output_dir mentioned
python $PWD/src/helpers/hyperparameters/find_chrombpnet_hyperparams.py \
       --genome=$reference_fasta \
       --bigwig=$bigwig_path \
       --peaks=$overlap_peak \
       --nonpeaks=$nonpeaks \
       --negative-sampling-ratio=$negative_sampling_ratio \
       --outlier_threshold=0.99 \
       --chr_fold_path=$fold \
       --inputlen=$inputlen \
       --outputlen=$outputlen \
       --filters=$filters \
       --n_dilation_layers=$n_dilation_layers \
       --bias_model_path=$bias_model \
       --output_dir=$output_dir 

# trains a bias model and store it in the directory mentioned with the given prefix
# make a check here that the model params are on same inout/output lenght
CUDA_VISIBLE_DEVICES=2 python $PWD/src/training/train.py \
       --genome=$reference_fasta \
       --bigwig=$bigwig_path \
       --peaks=$overlap_peak \
       --nonpeaks=$nonpeaks \
       --chr_fold_path=$fold \
       --params=$output_dir/chrombpnet_params.txt \
       --output_prefix=$output_dir/model.0 \
       --epochs=1 \
       --max_jitter=200 \
       --batch_size=64 \
       --architecture_from_file=$PWD/src/training/models/chrombpnet_model_with_bias_model.py \
       --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss 

# predictions and metrics on the bias model trained
CUDA_VISIBLE_DEVICES=2  python $PWD/src/training/predict.py \
        --genome=$reference_fasta \
        --bigwig=$bigwig_path \
        --peaks=$overlap_peak \
        --nonpeaks=$nonpeaks \
        --chr_fold_path=$fold \
        --inputlen=$inputlen \
        --outputlen=$outputlen \
        --output_prefix=$output_dir/ \
        --batch_size=64 \
        --model_h5=$output_dir/model.0.h5 \
