regions=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data_5M/peaks_no_blacklist.merged.bed
bigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_5M/GM12878_5M_07.18.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE/chrombpnet_model/interpret/full_GM12878_5M_wo_bias.bw
obigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_5M/GM12878_5M_07.18.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE/chrombpnet_model/interpret/wo_bias_expected_signal.txt

python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $bigwig --output_path $obigwig

head $obigwig

#head /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_5M/GM12878_5M_07.18.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE/chrombpnet_model/interpret/full_GM12878_5M.profile.stats

echo "\n"

regions=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data_5M/peaks_no_blacklist.merged.bed
bigwig=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data_5M/GM12878_unstranded.bw
obigwig=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data_5M/expected_signal_in_bw.txt

python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $bigwig --output_path $obigwig

head $obigwig
echo "\n"

