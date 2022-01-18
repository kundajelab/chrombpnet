dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/
CUDA_VISIBLE_DEVICES=1 python src/training/predict.py \
        --genome=$dir/reference/hg38.genome.fa \
        --bigwig=$dir/ATAC/ENCSR824DUE//preprocessing/bigWigs/ENCSR824DUE.bigWig \
        --peaks=$dir/ATAC/ENCSR824DUE//chrombpnet_model//filtered.peaks.bed \
        --nonpeaks=$dir/ATAC/ENCSR824DUE//chrombpnet_model//filtered.nonpeaks.bed \
        --chr_fold_path=$dir/splits/fold_0.json \
        --inputlen=2114 \
        --outputlen=1000 \
        --output_prefix=$dir/ATAC/ENCSR824DUE//chrombpnet_model//chrombpnet \
        --batch_size=256 \
        --model_h5=$dir/ATAC/ENCSR824DUE//chrombpnet_model//chrombpnet.h5

