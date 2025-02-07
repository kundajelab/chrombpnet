chrombpnet_nb=$1
chrombpnet=$2
cellline=$3
fold=$4
dtype=$5
gpu=$6


#regions=results/chrombpnet/auprc_curves/narrowpeak_genomewide_chr1.bed 
regions=results/chrombpnet/auprc_curves/downloads/$fold"_w_1000_s_250_narrowpeak.bed"
output_dir=results/chrombpnet/auprc_curves/$cellline/$dtype"_uncorrected"
mkdir $output_dir


chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
file=$output_dir/$cellline

echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig_no_bias.py  -cm $chrombpnet  --regions $regions -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline"_"$fold -t 1"
CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig_no_bias.py  -cm $chrombpnet  --regions $regions \
        -g $ref_fasta -b 32 -c $chrom_sizes -o $output_dir/$cellline"_"$fold -t 1

#echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig_new.py  -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline"_"$fold -t 1"
#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig_new.py  -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions \
#        -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline"_"$fold -t 1


#chrombpnet=results/chrombpnet/auprc_curves/$cellline/$dtype/$cellline"_"$fold"_w_bias_predictions.h5"
#chrombpnet_nb=results/chrombpnet/auprc_curves/$cellline/$dtype/$cellline"_"$fold"_wo_bias_predictions.h5"

#echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/make_only_bigwigs.py  -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline"_"$fold -t 1"
#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/make_only_bigwigs.py  -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions \
#        -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline"_"$fold -t 1



