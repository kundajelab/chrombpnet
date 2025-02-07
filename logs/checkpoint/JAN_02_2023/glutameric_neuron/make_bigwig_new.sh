chrombpnet_nb=$1
chrombpnet=$2
bias=$3
cellline=$4
gpu=$5
regions=$6
output_dir=$7

chrom_sizes=chrom.sizes
ref_fasta=hg38.genome.fa

#file=$output_dir/merged.$cellline
file=$output_dir/$cellline


#echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1"
#CUDA_VISIBLE_DEVICES=$gpu python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/predict_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions \
#	-g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1


#python normalize_bigwigs_from_peaks.py -bed merged.bed -bigwig $output_dir/$cellline"_wo_bias.bw" -o $output_dir/$cellline"_scale.txt"
#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/interpret/interpret.py -g $ref_fasta -r $regions -m $chrombpnet_nb -o $file --profile_or_counts counts
#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/interpret/interpret.py -g $ref_fasta -r $regions -m $bias -o $file --profile_or_counts counts

#python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.counts_scores.h5 -r $file.interpreted_regions_v2.bed -c $chrom_sizes -o $file.counts.bw -s $file.counts.stat -t 1

#python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.counts_scores.h5 -r $file.interpreted_regions.bed -c $chrom_sizes -o $file.counts.bw -s $file.counts.stat -t 1

#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/interpret/interpret.py -g $ref_fasta -r $regions -m $chrombpnet_nb -o $file --profile_or_counts profile
CUDA_VISIBLE_DEVICES=$gpu python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/interpret/interpret.py -g $ref_fasta -r $regions -m $bias -o $file --profile_or_counts counts
python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.counts_scores.h5 -r $file.interpreted_regions_v2.bed -c $chrom_sizes -o $file.counts.bw -s $file.counts.stat -t 1
CUDA_VISIBLE_DEVICES=$gpu python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/interpret/interpret.py -g $ref_fasta -r $regions -m $bias -o $file --profile_or_counts profile
python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.profile_scores.h5 -r $file.interpreted_regions_v2.bed -c $chrom_sizes -o $file.profile.bw -s $file.profile.stat -t 1

#python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.profile_scores.h5 -r $file.interpreted_regions.bed -c $chrom_sizes -o $file.profile.bw -s $file.profile.stat -t 1

#python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $output_dir/$cellline"bias.bw"  --output_path $file1"_scale_on_30k.txt"

