bigwig_path=$1
overlap_peak=$2
nonpeaks=$3
reference_fasta=$4
threshold_factor=$5
inputlen=$6
outputlen=$7
fold=$8

mkdir output
mkdir output/bias_model
output_dir=output/bias_model

#creates a bias_params.txt file in the output_dir mentioned
python $PWD/src/helpers/hyperparameters/find_bias_hyperparams.py \
       --bigwig=$bigwig_path \
       --peaks=$overlap_peak \
       --nonpeaks=$nonpeaks \
       --genome=$reference_fasta \
       --threshold_factor=$threshold_factor  \
       --chr_fold_path=$fold \
       --outputlen=$outputlen \
       --output_dir=$output_dir 

# trains a bias model and store it in the directory mentioned with the given prefix
CUDA_VISIBLE_DEVICES=0 python $PWD/src/training/train.py \
       --genome=$reference_fasta \
       --nonpeaks=$nonpeaks \
       --peaks=None \
       --output_prefix=$output_dir/model.0 \
       --chr_fold_path=$fold \
       --epochs=40 \
       --params=$output_dir/bias_params.txt \
       --inputlen=$inputlen \
       --outputlen=$outputlen \
       --max_jitter=10 \
       --batch_size=64 \
       --bigwig=$bigwig_path \
       --architecture_from_file=$PWD/src/training/models/bpnet_model.py \
       --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss 

# predictions and metrics on the bias model trained
CUDA_VISIBLE_DEVICES=0 python $PWD/src/training/predict.py \
        --genome=$reference_fasta \
        --peaks=$overlap_peak \
        --nonpeaks=$nonpeaks \
        --output_prefix=$output_dir \
        --chr_fold_path=$fold \
        --inputlen=$inputlen \
        --outputlen=$outputlen \
        --batch_size=64 \
        --model_h5=$output_dir/model.0.h5 \
        --params=$output_dir/bias_params.txt \
        --bigwig=$bigwig_path \
