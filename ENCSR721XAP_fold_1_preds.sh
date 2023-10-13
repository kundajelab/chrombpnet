chrombpnet_nb=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR721XAP/chrombppnet_model_encsr283tme_bias_fold_1/chrombpnet_wo_bias.h5
chrombpnet=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR721XAP/chrombppnet_model_encsr283tme_bias_fold_1/chrombpnet.h5
bias=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR721XAP/chrombppnet_model_encsr283tme_bias_fold_1/bias_model_scaled.h5
cellline=ENCSR721XAP
gpu=0


regions=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR721XAP/chrombppnet_model_encsr283tme_bias_fold_1/interpret/full_ENCSR721XAP.interpreted_regions_counts.bed
output_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/bigwigs/DNASE/ENCSR721XAP/chrombppnet_model_encsr283tme_bias_fold_1/preds/


chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
file=$output_dir/$cellline

echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1"
CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions \
        -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1
