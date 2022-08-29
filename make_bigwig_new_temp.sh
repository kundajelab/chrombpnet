chrombpnet_nb=$1
chrombpnet=$2
bias=$3
cellline=$4
gpu=$5
regions=$6
output_dir=$7

chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
#chrom_sizes=/mnt/data/annotations/by_release/hg19/hg19.chrom.sizes
#ref_fasta=/mnt/data/male.hg19.fa

file=$output_dir/merged.$cellline


#python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.counts_scores.h5 -r $file.interpreted_regions.bed -c $chrom_sizes -o $file.counts.bw -s $file.counts.stat -t 1


python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.profile_scores.h5 -r $file.interpreted_regions.bed -c $chrom_sizes -o $file.profile.bw -s $file.profile.stat -t 1
