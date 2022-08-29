reference_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
bigwig_path=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data/GM12878_unstranded.bw
#/oak/stanford/groups/akundaje/arpitas/thyroid_normal/croo/chrombpnet/data/downloads/merged-bigwig.bw
output_dir=$1
fold=$2
peaks_dir=$3
#/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/splits/fold_0.json
inputlen=2114
outputlen=1000
temp_dir=$output_dir/full_readdepth/
mkdir temp_dir
echo $output_dir/chrombpnet.h5
python $PWD/src/training/predict.py \
        --genome=$reference_fasta \
        --bigwig=$bigwig_path \
        --peaks=$peaks_dir/filtered.peaks.bed \
        --nonpeaks=$peaks_dir/filtered.nonpeaks.bed \
        --chr_fold_path=$fold \
        --inputlen=$inputlen \
        --outputlen=$outputlen \
        --output_prefix=$temp_dir/chrombpnet \
        --batch_size=256 \
        --model_h5=$output_dir/chrombpnet.h5 | tee -a $logfile

