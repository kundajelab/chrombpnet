chrombpnet_nb=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/ATAC/ENCSR483RKN/chrombpnet_model_feb15/chrombpnet_wo_bias.h5
chrombpnet=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/ATAC/ENCSR483RKN/chrombpnet_model_feb15/chrombpnet.h5
bias=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/ATAC/ENCSR483RKN/chrombpnet_model_feb15/bias_model_scaled.h5
cellline=K562
gpu=0

regions=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/ATAC/ENCSR483RKN/negatives_data/negatives_with_summit_subsample.bed
output_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/bigwigs/ATAC/ENCSR483RKN/background/
mkdir $output_dir


chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
file=$output_dir/$cellline

echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1"
CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions \
	-g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1


