dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/bias_model/interpret/
file=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/benchmarking/bias/GM12878.bias
chrom_sizes=reference/chrom.sizes

#python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $dir/GM12878.profile_scores.h5 -r $dir/GM12878.interpreted_regions.bed -c $chrom_sizes -o $file.profile.bw -s $file.profile.stat -t 1

regions=$dir/GM12878.interpreted_regions.bed
bigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/benchmarking/chrombpnet/chrombpnetbias.bw
obigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/benchmarking/bias/GM12878_bias_scaled.bw


#python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $bigwig --output_path $obigwig
#python normalize_shap_bigwig_tracks.py --input_stat $file.profile.stat -bigwig $file.profile.bw -o $file.profile.scaled.bw

regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/uncorrected_model_05.10.2022/uncorrected_model/interpret/GM12878.interpreted_regions.merged.bed
python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $bigwig --output_path $obigwig


