chrombpnet_nb=$1
cellline=$2
model_dir=$3
gpu=$4
ddtype=$5 

regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$ddtype/$cellline/negatives_data/negatives_subsample.bed
output_dir=$model_dir/background_interpret_new/
mkdir $output_dir


chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
file=$output_dir/$cellline



CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/interpret/interpret.py -g $ref_fasta -r $regions -m $chrombpnet_nb -o $file --profile_or_counts counts
python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.counts_scores.h5 -r $file.interpreted_regions_v2.bed -c $chrom_sizes -o $file.counts.bw -s $file.counts.stat -t 1

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/interpret/interpret.py -g $ref_fasta -r $regions -m $chrombpnet_nb -o $file --profile_or_counts profile
python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.profile_scores.h5 -r $file.interpreted_regions_v2.bed -c $chrom_sizes -o $file.profile.bw -s $file.profile.stat -t 1


