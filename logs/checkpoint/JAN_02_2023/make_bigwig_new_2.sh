chrombpnet_nb=$1
chrombpnet=$2
bias=$3
cellline=$4
gpu=$5
regions=$6
output_dir=$7
file1=$8

chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
#chrom_sizes=/mnt/data/annotations/by_release/hg19/hg19.chrom.sizes
#ref_fasta=/mnt/data/male.hg19.fa

file=$output_dir/merged.$cellline


regions=$file1.interpreted_regions.bed
#echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1"
#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions \
#	-g $ref_fasta -c $chrom_sizes -o $file1 -t 1


#python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file1.counts_scores.h5 -r $regions -c $chrom_sizes -o $file1.counts.bw -s $file1.counts.stat -t 1


#python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file1.profile_scores.h5 -r $regions -c $chrom_sizes -o $file1.profile.bw -s $file1.profile.stat -t 1

regions=$file1.interpreted_regions_v2.merged.bed

#python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $file1"_wo_bias.bw" --output_path $file1"_wo_bias_scaled.txt"
python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig results/chrombpnet/ATAC_PE/GM12878/data/GM12878_unstranded.bw  --output_path $file1"_scale_on_30k.txt"
#python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $file1"bias.bw" --output_path $file1"_bias_scaled.txt"
#python normalize_shap_bigwig_tracks.py --input_stat $output_dir/merged.$cellline.profile.stat -bigwig $output_dir/merged.$cellline.profile.bw -o $output_dir/merged.$cellline.profile.scaled.bw
