reference_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
bigwig_path=/oak/stanford/groups/akundaje/arpitas/thyroid_normal/croo/chrombpnet/data/downloads/merged-bigwig.bw
output_dir=/oak/stanford/groups/akundaje/arpitas/thyroid_normal/croo/chrombpnet/models/chrombpnet_default_bias_model/
fold=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/splits/fold_0.json
inputlen=2114
outputlen=1000
temp_dir=temp
echo $output_dir/chrombpnet.h5
python $PWD/src/training/predict.py \
        --genome=$reference_fasta \
        --bigwig=$bigwig_path \
        --peaks=$output_dir/filtered.peaks.bed \
        --nonpeaks=$output_dir/filtered.nonpeaks.bed \
        --chr_fold_path=$fold \
        --inputlen=$inputlen \
        --outputlen=$outputlen \
        --output_prefix=$temp_dir/chrombpnet \
        --batch_size=256 \
        --model_h5=$output_dir/chrombpnet.h5 | tee -a $logfile

