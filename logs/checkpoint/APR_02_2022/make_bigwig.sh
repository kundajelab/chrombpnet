chrombpnet_nb=$1
chrombpnet=$2
bias=$3
cellline=$4
gpu=$5
regions=$6
output_dir=$7

chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa


file=$output_dir/$cellline

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/interpret/interpret.py -g $ref_fasta -r $regions -m $chrombpnet_nb -o $file

echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1"
CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions \
	-g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1



python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.counts_scores.h5 -r $output_dir/$cellline.interpreted_regions.bed -c $chrom_sizes -o $file.counts.bw -s $file.counts.stat -t 1
python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.profile_scores.h5 -r $output_dir/$cellline.interpreted_regions.bed -c $chrom_sizes -o $file.profile.bw -s $file.profile.stat -t 1

