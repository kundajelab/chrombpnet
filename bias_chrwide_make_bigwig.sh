chrombpnet_nb=$1
chrombpnet=$2
cellline=$3
fold=$4
dtype=$5
gpu=$6


#regions=results/chrombpnet/auprc_curves/narrowpeak_genomewide_chr1.bed 
regions=results/chrombpnet/auprc_curves/downloads/$fold"_w_1000_s_250_narrowpeak.bed"
output_dir=/mnt/lab_data2/anusri/print_analysis/make_bg_bigwig/$cellline/
#=results/chrombpnet/auprc_curves/$cellline/$dtype"_uncorrected"
mkdir $output_dir

echo "$output_dir"
chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
file=$output_dir/$cellline

#predict_to_bigwig_print_bias.py

regions=results/chrombpnet/auprc_curves/downloads/$fold"_w_1000_s_250_narrowpeak.bed"

echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig_no_bias.py  -cm $chrombpnet  --regions $regions -g $ref_fasta -c $chrom_sizes -o  $output_dir/"ATAC_bias_"$fold -t 1"
CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig_no_bias.py  -cm $chrombpnet  --regions $regions \
        -g $ref_fasta -b 32 -c $chrom_sizes -o $output_dir/"ATAC_bias_"$fold -t 1

fold=fold1
regions=results/chrombpnet/auprc_curves/downloads/$fold"_w_1000_s_250_narrowpeak.bed"

echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig_no_bias.py  -cm $chrombpnet  --regions $regions -g $ref_fasta -c $chrom_sizes -o  $output_dir/"ATAC_bias_"$fold -t 1"
CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig_no_bias.py  -cm $chrombpnet  --regions $regions \
        -g $ref_fasta -b 32 -c $chrom_sizes -o $output_dir/"ATAC_bias_"$fold -t 1


fold=fold2
regions=results/chrombpnet/auprc_curves/downloads/$fold"_w_1000_s_250_narrowpeak.bed"

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig_no_bias.py  -cm $chrombpnet  --regions $regions \
        -g $ref_fasta -b 32 -c $chrom_sizes -o $output_dir/"ATAC_bias_"$fold -t 1


fold=fold3
regions=results/chrombpnet/auprc_curves/downloads/$fold"_w_1000_s_250_narrowpeak.bed"

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig_no_bias.py  -cm $chrombpnet  --regions $regions \
        -g $ref_fasta -b 32 -c $chrom_sizes -o $output_dir/"ATAC_bias_"$fold -t 1

fold=fold4
regions=results/chrombpnet/auprc_curves/downloads/$fold"_w_1000_s_250_narrowpeak.bed"

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig_no_bias.py  -cm $chrombpnet  --regions $regions \
        -g $ref_fasta -b 32 -c $chrom_sizes -o $output_dir/"ATAC_bias_"$fold -t 1

