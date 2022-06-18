#regions=results/chrombpnet/ATAC_PE/K562/ATAC_PE_12.30.2021/chrombpnet_model/testing/beta_globin.coords
regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/data/chr20.peaks.bed
#regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/data/chr20.peaks.bed
#regions=/srv/scratch/anusri/colloboration_data_overlap/interesting.locus


output_dir=results/chrombpnet/ATAC_PE/K562/ATAC_PE_01.24.2022_bigger_jitter_on_bias/chrombpnet_model/chr20_interpret/
moutput_dir=results/chrombpnet/ATAC_PE/K562/ATAC_PE_01.24.2022_bigger_jitter_on_bias/chrombpnet_model/
mkdir $output_dir
gpu=0


python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
        -r $regions \
	-h5 $output_dir/chrombpnet_wo_bias.counts_scores.h5 \
	-o $output_dir/chrombpnet_wo_bias_counts.bw \
	-s $output_dir/chrombpnet_wo_bias_counts.stats \
	-t 1


