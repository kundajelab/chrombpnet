model_h5=results/chrombpnet/ATAC_PE/IMR90/nautilus_runs_apr12/IMR90_04.09.2022_bias_128_4_1234_0.4_fold_0/bias_model/bias.h5
chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes
reference_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
fold=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/splits/fold_4.json
bigwig_path=results/chrombpnet/ATAC_PE/IMR90/data/IMR90_unstranded.bw


CUDA_VISIBLE_DEVICES=5 python $PWD/src/training/predict.py \
        --genome=$reference_fasta \
        --bigwig=$bigwig_path \
        --nonpeaks=results/chrombpnet/ATAC_PE/IMR90/nautilus_runs_apr12/IMR90_04.09.2022_bias_128_4_1234_0.4_fold_0/bias_model/filtered.bias_nonpeaks.bed \
        --chr_fold_path=$fold \
        --inputlen=2114 \
        --outputlen=1000 \
        --output_prefix=bias_model_4/ \
        --batch_size=256 \
        --model_h5=$model_h5 \



