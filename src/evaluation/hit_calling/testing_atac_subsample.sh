output_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/$1/$2/09_06_2022_motif_scanning/moods_run/
peak_bed=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/interpret/merged.GM12878.interpreted_regions.bed
modisco_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/GM12878/modisco_crop_500_100K_seqs_1
shap_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/$1/$2/chrombpnet_model/interpret/full_$1
genome=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa 

mkdir  $output_dir
echo $output_dir

# update to the new peak file
peak_bed=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/09_06_2022_motif_scanning/mooods_run/equal.spaced.merged.k562.bed


bash run_moods_hits_test_subsample.sh $modisco_dir $shap_dir $output_dir $peak_bed "mean_norm"

