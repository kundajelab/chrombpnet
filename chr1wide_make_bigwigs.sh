chrombpnet_nb=$1
chrombpnet=$2
cellline=$3
gpu=$4

regions=results/chrombpnet/auprc_curves/narrowpeak_genomewide_chr1.bed 
output_dir=results/chrombpnet/auprc_curves/$cellline
mkdir $output_dir


chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
file=$output_dir/$cellline

echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig_new.py  -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1"
CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig_new.py  -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions \
        -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1


chrombpnet=results/chrombpnet/auprc_curves/$cellline/$cellline"_w_bias_predictions.h5"
chrombpnet_nb=results/chrombpnet/auprc_curves/$cellline/$cellline"_wo_bias_predictions.h5"

echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/make_only_bigwigs.py  -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1"
CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/make_only_bigwigs.py  -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions \
        -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1



