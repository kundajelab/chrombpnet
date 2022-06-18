output_dir=results/chrombpnet/ATAC_PE/K562/ATAC_PE_12.30.2021/chrombpnet_model/interpret
regions=$output_dir/K562.interpreted_regions.bed
cell_type=K562

output_dir_n=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/ATAC_PE/K562/ATAC_PE_12.30.2021/SIGNAL/

python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
        -r $regions \
	-h5 $output_dir_n/K562.profile_scores.h5 \
	-o $output_dir_n/K562_profile.bw \
	-s $output_dir_n/K562_profile.stats \
	-t 1

python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
        -r $regions \
	-h5 $output_dir_n/K562.counts_scores.h5 \
	-o $output_dir_n/K562_counts.bw \
	-s $output_dir_n/K562_counts.stats \
	-t 1



