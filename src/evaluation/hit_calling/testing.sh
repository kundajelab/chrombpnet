bed_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/nautilus_runs_may18/K562_05.13.2022_bias_128_4_1234_0.5_fold_0/interpret/merged.K562.interpreted_regions.bed

#bedtools slop -i /mnt/data/annotations/blacklist/GRch38/GRch38_unified_blacklist.bed.gz -g /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes -b 1057 > temp.bed
#bedtools intersect -v -a $bed_dir -b temp.bed | awk ' {print $1 "\t" $2 + $10 - 500 "\t" $2 + $10 + 500}' > merged_peaks_no_lacklist.w1000.bed
genom=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa 
 
#bedtools sort -i merged_peaks_no_lacklist.w1000.bed | bedtools merge -i stdin > merged.k562.dnase.bed

#python create_equal_width_peaks.py --bed merged.k562.dnase.bed --outf equal.spaced.merged.k562.new.bed

peaks_bed=equal.spaced.merged.k562.new.bed

wc -l equal.spaced.merged.k562.new.bed

peak_bed=equal.spaced.merged.k562.new.bed

modisco_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/DNASE/K562/modisco_crop_500_100K_seqs_1
shap_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/nautilus_runs_may18/K562_05.13.2022_bias_128_4_1234_0.5_fold_0/interpret/merged.K562
#output_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/K562_02.08.2022_bias_128_4_2356/06_22_2022_motif_scanning/mean_moods_baseairmodels_counts_new_divide_by_total_score_pval_0.01_moods/
output_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/K562_02.08.2022_bias_128_4_2356/06_22_2022_motif_scanning/mean_moods_baseairmodels_counts_new_divide_by_total_score/
mkdir  $output_dir
echo $output_dir

#bash run_moods_hits.sh $modisco_dir $shap_dir $output_dir $peak_bed "mean_norm"
#python break_overlaps.py -o $output_dir
python break_overlaps_with_cwm_score_2.py -o $output_dir -tfm $modisco_dir//modisco_results_allChroms_counts.hdf5 -bw $shap_dir.counts.bw -g $genom




modisco_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/DNASE/K562/modisco_crop_500_100K_seqs_1
shap_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/nautilus_runs_may18/K562_05.13.2022_bias_128_4_1234_0.5_fold_0/interpret/merged.K562
output_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/K562_02.08.2022_bias_128_4_2356/06_22_2022_motif_scanning/mean_moods_baseairmodels_counts_new/
mkdir  $output_dir
echo $output_dir

#bash run_moods_hits.sh $modisco_dir $shap_dir $output_dir $peak_bed "mean_sum"
#python break_overlaps.py -o $output_dir
#python break_overlaps_with_cwm_score_2.py -o $output_dir -tfm $modisco_dir//modisco_results_allChroms_counts.hdf5 -bw $shap_dir.counts.bw -g $genom




modisco_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/DNASE/K562/modisco_crop_500_100K_seqs_1
shap_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/nautilus_runs_may18/K562_05.13.2022_bias_128_4_1234_0.5_fold_0/interpret/merged.K562
output_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/K562_02.08.2022_bias_128_4_2356/06_22_2022_motif_scanning/moods_baseairmodels_counts_new_divide_by_total_score/
mkdir  $output_dir
echo $output_dir

#bash run_moods_hits.sh $modisco_dir $shap_dir $output_dir $peak_bed "norm"
#python break_overlaps.py -o $output_dir
#python break_overlaps_with_cwm_score_2.py -o $output_dir -tfm $modisco_dir//modisco_results_allChroms_counts.hdf5 -bw $shap_dir.counts.bw -g $genom


modisco_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/DNASE/K562/modisco_crop_500_100K_seqs_1
shap_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/nautilus_runs_may18/K562_05.13.2022_bias_128_4_1234_0.5_fold_0/interpret/merged.K562
output_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/K562_02.08.2022_bias_128_4_2356/06_22_2022_motif_scanning/moods_baseairmodels_counts_new/
mkdir  $output_dir
echo $output_dir

#bash run_moods_hits.sh $modisco_dir $shap_dir $output_dir $peak_bed "sum"
#python break_overlaps.py -o $output_dir
#python break_overlaps_with_cwm_score_2.py -o $output_dir -tfm $modisco_dir//modisco_results_allChroms_counts.hdf5 -bw $shap_dir.counts.bw -g $genom


